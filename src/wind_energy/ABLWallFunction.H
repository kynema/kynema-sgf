#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/wind_energy/MOData.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/** Wall-function computations for LES simulations of ABL
 *  \ingroup we_abl
 *
 *  This class performs the necessary computations at the beginning of
 *  predictor/corrector steps. The actual BC population in ghost cells is
 *  performed by ABLVelWallFunc BC interface class.
 */
class ABLWallFunction
{
public:
    explicit ABLWallFunction(const CFDSim& sim);

    ~ABLWallFunction() = default;

    MOData& mo() { return m_mo; }
    [[nodiscard]] const MOData& mo() const { return m_mo; }

    //! Return the plane-averaged computed friction velocity at any given
    //! instance
    [[nodiscard]] amrex::Real utau() const { return m_mo.utau; }

    //! Initialize the log-law height based on user inputs
    void init_log_law_height();

    //! Update the mean velocity at a given timestep
    void update_umean(
        const VelPlaneAveragingFine& vpa, const FieldPlaneAveragingFine& tpa);

    void update_tflux(amrex::Real tflux);

private:
    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    //! Monin-Obukhov instance
    MOData m_mo;

    int m_direction{2};   ///< Direction normal to wall
    bool m_use_fch{true}; ///< Use first cell height?

    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    //! Ability to read in a table of surface temperature versus time.
    std::string m_surf_temp_timetable;
    amrex::Vector<amrex::Real> m_surf_temp_time;
    amrex::Vector<amrex::Real> m_surf_temp_value;

    bool m_tempflux{true};
    bool m_temp_table{false};
    amrex::Real m_surf_temp_rate{0.0_rt};
    amrex::Real m_surf_temp_rate_tstart{0.0_rt};
    amrex::Real m_surf_temp_init{300.0_rt};

    bool m_inflow_outflow{false};
    amrex::Real m_wf_vmag{0.0_rt};
    amrex::Array<amrex::Real, 2> m_wf_vel{0.0_rt, 0.0_rt};
    amrex::Real m_wf_theta{300.0_rt};
};

/** Applies a shear-stress value at the domain boundary
 *  \ingroup field_bc we_abl
 *
 *  \sa ABLWallFunction
 */
class ABLVelWallFunc : public FieldBCIface
{
public:
    ABLVelWallFunc(Field& velocity, const ABLWallFunction& wall_func);

    void operator()(Field& velocity, FieldState rho_state) override;

    template <typename ShearStress>
    void
    wall_model(Field& velocity, FieldState rho_state, const ShearStress& tau);

private:
    const ABLWallFunction& m_wall_func;
    std::string m_wall_shear_stress_type{"moeng"};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
};

class ABLTempWallFunc : public FieldBCIface
{
public:
    ABLTempWallFunc(Field& temperature, const ABLWallFunction& wall_fuc);

    void operator()(Field& temperature, FieldState rho_state) override;

    template <typename HeatFlux>
    void
    wall_model(Field& temperature, FieldState rho_state, const HeatFlux& tau);

private:
    const ABLWallFunction& m_wall_func;
    std::string m_wall_shear_stress_type{"moeng"};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
};

} // namespace amr_wind

#endif /* ABLWALLFUNCTION_H */
