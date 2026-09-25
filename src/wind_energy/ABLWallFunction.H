#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "src/CFDSim.H"
#include "src/utilities/FieldPlaneAveraging.H"
#include "src/core/FieldBCOps.H"
#include "src/wind_energy/MOData.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

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

    //! Monin-Obukhov data carrying the mean quantities of a given mesh level
    //!
    //! On a single-level mesh this is the same object as mo(). When more
    //! than one level touches the wall, every level holds the plane averages
    //! of its own wall-adjacent cells at their own height, so that the wall
    //! models that combine the local first-cell values with their plane
    //! averages return the same mean stress and heat flux on every level.
    //! The friction velocity, Obukhov length and surface heat flux are
    //! shared by all levels.
    //!
    //! \param lev [in] Mesh level
    [[nodiscard]] const MOData& mo(int lev) const
    {
        return (lev < static_cast<int>(m_mo_lev.size())) ? m_mo_lev[lev] : m_mo;
    }

    //! Return the plane-averaged computed friction velocity at any given
    //! instance
    [[nodiscard]] amrex::Real utau() const { return m_mo.utau; }

    //! Initialize the log-law height based on user inputs
    void init_log_law_height(int max_level);

    //! Update the mean velocity at a given timestep
    void
    update_umean(const VelPlaneAveraging& vpa, const FieldPlaneAveraging& tpa);

    void update_tflux(amrex::Real tflux);

    //! Compute the mean quantities of the wall-adjacent cells of every level
    //! that touches the wall (multi-level meshes only).
    //! Public because it launches GPU kernels
    void update_level_means();

private:
    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    //! Monin-Obukhov instance
    MOData m_mo;

    //! Monin-Obukhov data per mesh level, empty on single-level meshes
    amrex::Vector<MOData> m_mo_lev;

    int m_direction{2};             ///< Direction normal to wall
    bool m_use_fch{true};           ///< Use first cell height?
    amrex::Real m_wall_pos{0.0_rt}; ///< Position of the wall
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

    //! Apply the shear-stress model, built per level from the mean
    //! quantities of that level
    template <typename ShearStress>
    void wall_model(Field& velocity, FieldState rho_state);

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

    //! Apply the heat-flux model, built per level from the mean quantities
    //! of that level
    template <typename HeatFlux>
    void wall_model(Field& temperature, FieldState rho_state);

private:
    const ABLWallFunction& m_wall_func;
    std::string m_wall_shear_stress_type{"moeng"};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
};

} // namespace kynema_sgf

#endif /* ABLWALLFUNCTION_H */
