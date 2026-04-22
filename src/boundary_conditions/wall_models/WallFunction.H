#ifndef WALLFUNCTION_H
#define WALLFUNCTION_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "amr-wind/boundary_conditions/wall_models/LogLaw.H"
#include "amr-wind/boundary_conditions/wall_models/MOSD.H"

namespace amr_wind {

/** Wall-function computations for LES simulations
 *  \ingroup utilities
 *
 *  This class performs the necessary computations at the beginning of
 *  predictor/corrector steps. The actual BC population in ghost cells is
 *  performed by VelWallFunc BC interface class.
 */
class WallFunction
{
public:
    explicit WallFunction(CFDSim& sim);

    [[nodiscard]] amrex::Real utau() const { return m_log_law.utau_mean; }
    [[nodiscard]] LogLaw log_law() const { return m_log_law; }
    [[nodiscard]] MOSD mosd() const { return m_mosd; }

    //! Update the mean velocity at a given timestep
    void update_umean();
    void update_utau_mean();
    void update_time();

    ~WallFunction() = default;

private:
    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    //! LogLaw instance
    LogLaw m_log_law;
    int m_direction{2}; ///< Direction normal to wall, hardcoded to z
    VelPlaneAveragingFine m_pa_vel;

    MOSD m_mosd;
};

/** Applies a shear-stress value at the domain boundary
 *  \ingroup field_bc utilities
 *
 *  \sa WallFunction
 */
class VelWallFunc : public FieldBCIface
{
public:
    VelWallFunc(Field& velocity, WallFunction& wall_func);

    void operator()(Field& velocity, FieldState rho_state) override;

    static void
    wall_model(Field& velocity, FieldState rho_state, amrex::Real utau);

    template <typename ShearStressSimple>
    static void wall_model(
        Field& velocity, FieldState rho_state, const ShearStressSimple& tau);

private:
    WallFunction& m_wall_func;
    std::string m_wall_shear_stress_type{"constant"};
};
} // namespace amr_wind

#endif /* WALLFUNCTION_H */
