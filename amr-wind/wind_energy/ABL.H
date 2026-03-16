#ifndef ABL_H
#define ABL_H

#include <memory>

#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/wind_energy/ABLFieldInitFile.H"
#include "amr-wind/wind_energy/ABLWallFunction.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/wind_energy/ABLStatsBase.H"
#include "amr-wind/wind_energy/ABLModulatedPowerLaw.H"
#include "amr-wind/wind_energy/ABLAnelastic.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/**
 *  \defgroup wind Wind energy modeling
 *  Wind energy modeling
 *
 *  This group documents all the wind-energy related physics models available in
 *  AMR-Wind
 *
 *  \ingroup physics
 */

/**
 *  \defgroup we_abl ABL
 *  Atmospheric boundary layer (ABL) modeling
 *
 *  \ingroup wind
 */

class incflo;

namespace amr_wind {

class ABLMesoscaleInput;

namespace pde::icns {
class ABLForcing;
class ABLMeanBoussinesq;
class ABLMesoForcingMom;
class HurricaneForcing;
} // namespace pde::icns

namespace pde::temperature {
class ABLMesoForcingTemp;
class HurricaneTempForcing;
} // namespace pde::temperature

/** Atmospheric Boundary Layer physics
 *  \ingroup we_abl
 *
 *  This class implements the necessary code to simulation atmospheric boundary
 *  layers with AMR-Wind. It reads the namespace `ABL` from the input file to
 *  determine how the initial conditions are generated. At each timestep it
 *  performs plane averages to collect statistics as well as set the necessary
 *  data to be used by wall models as well as any ABLForcing terms that might be
 *  active.
 *
 *  \sa ABLFieldInit, ABLWallFunction, FieldPlaneAveraging, ABLBoundaryPlane
 */
class ABL : public Physics::Register<ABL>
{
public:
    static std::string identifier() { return "ABL"; }

    explicit ABL(CFDSim& sim);

    ~ABL() override;

    const ABLWallFunction& abl_wall_function() const { return m_abl_wall_func; }

    void post_init_actions() override;

    void post_regrid_actions() override;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_advance_work() override;

    void pre_predictor_work() override;

    void post_advance_work() override;

    void register_forcing_term(pde::icns::ABLForcing* forcing) const
    {
        m_abl_forcing = forcing;
    }

    void register_mean_boussinesq_term(pde::icns::ABLMeanBoussinesq* term) const
    {
        m_abl_mean_bous = term;
    }
    void register_meso_mom_forcing(pde::icns::ABLMesoForcingMom* forcing) const
    {
        m_abl_meso_mom_forcing = forcing;
    }
    void register_meso_temp_forcing(
        pde::temperature::ABLMesoForcingTemp* forcing) const
    {
        m_abl_meso_temp_forcing = forcing;
    }

    void register_hurricane_forcing(pde::icns::HurricaneForcing* forcing) const
    {
        m_hurricane_forcing = forcing;
    }

    void register_hurricane_temp_forcing(
        pde::temperature::HurricaneTempForcing* forcing) const
    {
        m_hurricane_temp_forcing = forcing;
    }

    const ABLBoundaryPlane& bndry_plane() const { return *m_bndry_plane; }
    ABLBoundaryPlane& bndry_plane() { return *m_bndry_plane; }
    const ABLModulatedPowerLaw& abl_mpl() const { return *m_abl_mpl; }

    //! Return the ABL statistics calculator
    const ABLStatsBase& abl_statistics() const { return *m_stats; }

    const ABLMesoscaleInput& abl_meso_file() const { return *m_meso_file; }

    const ABLAnelastic& anelastic() const { return *m_abl_anelastic; }

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_mueff;
    Field& m_density;
    Field* m_temperature{nullptr};
    Field* m_tke{nullptr};
    Field* m_sdr{nullptr};

    ABLWallFunction m_abl_wall_func;

    //! ABL field initializer instance
    std::unique_ptr<ABLFieldInit> m_field_init;

    //! ABL boundary plane instance
    std::unique_ptr<ABLBoundaryPlane> m_bndry_plane;

    //! ABL modulated power law instance
    std::unique_ptr<ABLModulatedPowerLaw> m_abl_mpl;

    //! Reference to ABL forcing term if present
    mutable pde::icns::ABLForcing* m_abl_forcing{nullptr};

    //! ABL integrated statistics object
    std::unique_ptr<ABLStatsBase> m_stats;

    //! ABL anelastic instance
    std::unique_ptr<ABLAnelastic> m_abl_anelastic;

    mutable pde::icns::ABLMeanBoussinesq* m_abl_mean_bous{nullptr};

    mutable pde::icns::ABLMesoForcingMom* m_abl_meso_mom_forcing{nullptr};
    mutable pde::temperature::ABLMesoForcingTemp* m_abl_meso_temp_forcing{
        nullptr};

    std::unique_ptr<ABLMesoscaleInput> m_meso_file;

    mutable pde::icns::HurricaneForcing* m_hurricane_forcing{nullptr};
    mutable pde::temperature::HurricaneTempForcing* m_hurricane_temp_forcing{
        nullptr};
    //! Default value set based on https://turbmodels.larc.nasa.gov/sst.html
    amrex::Real m_init_sdr{25.0_rt};

    //! Hybrid RANS-LES with Nalu-wind
    bool m_hybrid_rl{false};

    //! File input flag
    bool m_file_input{false};
    //! File input initializer
    std::unique_ptr<ABLFieldInitFile> m_field_init_file;
};

} // namespace amr_wind

#endif /* ABL_H */
