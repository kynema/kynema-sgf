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
    //! On a single-level mesh, when only level 0 touches the wall, when
    //! the bottom boundary is not a wall model, in inflow-outflow mode, or with
    //! ABL.level_means_method = none, this is the same object as mo(). When a
    //! finer level touches the wall, every level with wall cells whose first
    //! cell sits above the roughness height holds the plane averages of the
    //! reference height carried to its own first-cell height along the
    //! Monin-Obukhov profile (any other level holds a copy of mo()), so that
    //! the wall models that combine the local first-cell values with their
    //! plane averages return the same mean stress and heat flux on every level
    //! when the plane-mean profiles follow that profile (see
    //! ABL.level_means_method). The friction velocity, Obukhov length and
    //! surface heat flux are shared by all levels.
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

private:
    //! Set the Monin-Obukhov data of every level that touches the wall (only
    //! when a level finer than level 0 does), see mo(int)
    void update_level_means();

    const CFDSim& m_sim;

    const amrex::AmrCore& m_mesh;

    //! Monin-Obukhov instance
    MOData m_mo;

    //! Monin-Obukhov data per mesh level; empty when mo(lev) is mo() on
    //! every level (see mo(int))
    amrex::Vector<MOData> m_mo_lev;

    //! How the per-level data are obtained: "mo_profile" (reference-height
    //! plane averages carried along the Monin-Obukhov profile) or "none"
    //! (every level uses the reference-height plane averages)
    std::string m_level_means_method{"mo_profile"};

    //! Whether the warning about a first cell at or below the roughness
    //! height was printed
    bool m_warned_low_first_cell{false};

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
    //! quantities of that level (from those of the reference height for the
    //! Donelan model)
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
    //! of that level (from those of the reference height for the Donelan
    //! model)
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
