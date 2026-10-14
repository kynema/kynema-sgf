#ifndef ABLFORCING_H
#define ABLFORCING_H

#include <limits>

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/utilities/trig_ops.H"
#include "src/utilities/linear_interpolation.H"
#include "src/utilities/FieldPlaneAveraging.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Forcing term to drive the ABL profile to a desired state
 *  \ingroup icns_src we_abl
 *
 *  The forcing is a uniform horizontal acceleration \f$\mathbf{F}\f$ that
 *  drives the planar-averaged velocity at `abl_forcing_height` to the target
 *  velocity. With `free_atmosphere_damping` enabled, the planar-averaged
 *  velocity above the free-atmosphere height \f$h\f$ is also relaxed toward
 *  the geostrophic wind that balances \f$\mathbf{F}\f$,
 *
 *  \f[
 *    U_g = F_y / f, \quad V_g = -F_x / f, \quad
 *    S_{x,y} = \frac{1}{\tau}
 *      \left( U_{g,x,y} - \langle u_{x,y} \rangle (z) \right)
 *      \quad \text{for} \quad z \geq h,
 *  \f]
 *
 *  where \f$f = 2 \Omega \sin \phi\f$ comes from the CoriolisForcing inputs.
 *  This removes the inertial oscillations the free atmosphere undergoes when
 *  the forcing changes in time.
 *
 *  \sa ABL
 */
class ABLForcing : public MomentumSource::Register<ABLForcing>
{
public:
    static std::string identifier() { return "ABLForcing"; }

    explicit ABLForcing(const CFDSim& sim);

    ~ABLForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

    void set_target_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_target_vel[0] = ux;
        m_target_vel[1] = uy;
    }

    void set_mean_velocities(amrex::Real ux, amrex::Real uy)
    {
        m_mean_vel[0] = ux;
        m_mean_vel[1] = uy;

        const auto& current_time = m_time.current_time();
        const auto& new_time = m_time.new_time();
        const auto& nph_time = 0.5_rt * (current_time + new_time);
        const auto& dt = m_time.delta_t();
        const auto& t_step = m_time.time_index();

        if (!m_vel_timetable.empty()) {
            // Forces should be applied at n+1/2. Because ABL forcing is a
            // delta, the difference between the target velocity (at n+1) and
            // the current velocity (at n) puts the force term at n+1/2
            const amrex::Real new_spd = ::kynema_sgf::interp::linear(
                m_time_table, m_speed_table, new_time);
            const amrex::Real new_dir = ::kynema_sgf::interp::linear_angle(
                m_time_table, m_direction_table, new_time,
                2.0_rt * std::numbers::pi_v<amrex::Real>);

            m_target_vel[0] = new_spd * std::cos(new_dir);
            m_target_vel[1] = new_spd * std::sin(new_dir);
        }

        m_abl_forcing[0] = (m_target_vel[0] - m_mean_vel[0]) / dt;
        m_abl_forcing[1] = (m_target_vel[1] - m_mean_vel[1]) / dt;

        if (m_write_force_timetable &&
            amrex::ParallelDescriptor::IOProcessor() &&
            (t_step % m_force_outfreq == 0) &&
            (current_time >= m_force_outstart)) {
            std::ofstream outfile;
            // Forces are recorded at n+1/2
            outfile.open(m_force_timetable, std::ios::out | std::ios::app);
            outfile << std::setprecision(17) << nph_time << "\t"
                    << m_abl_forcing[0] << "\t" << m_abl_forcing[1] << "\t"
                    << 0.0_rt;
            if (m_fa_damping) {
                const auto ug = geostrophic_velocity();
                outfile << "\t" << ug[0] << "\t" << ug[1] << "\t"
                        << m_fa_height;
            }
            outfile << '\n';
        }
    }

    [[nodiscard]] amrex::RealArray abl_forcing() const { return m_abl_forcing; }

    //! True if the free atmosphere is damped toward the geostrophic wind
    [[nodiscard]] bool free_atmosphere_damping() const { return m_fa_damping; }

    //! True if the free-atmosphere height follows the capping inversion
    [[nodiscard]] bool detect_free_atmosphere_height() const
    {
        return m_fa_detect_height;
    }

    //! Height above which the free-atmosphere damping acts
    [[nodiscard]] amrex::Real free_atmosphere_height() const
    {
        return m_fa_height;
    }

    //! Geostrophic wind that balances the current forcing
    [[nodiscard]] amrex::Array<amrex::Real, 2> geostrophic_velocity() const
    {
        return {
            m_abl_forcing[1] / m_fa_coriolis_factor,
            -m_abl_forcing[0] / m_fa_coriolis_factor};
    }

    /** Update the free-atmosphere damping inputs at the start of a timestep
     *
     *  The damping acts from the second update on: the first forcing after
     *  start-up corrects no drift yet, so it does not give the geostrophic
     *  wind.
     *
     *  \param vpa Planar-averaged velocity profile
     *  \param zi Capping inversion height above the bottom of the domain,
     *  used only when `detect_free_atmosphere_height` is on
     */
    void update_free_atmosphere(const VelPlaneAveraging& vpa, amrex::Real zi);

    /** Add the free-atmosphere damping to the source term
     *
     *  Public because it launches a device lambda.
     *
     *  \param lev Mesh level
     *  \param src_term Momentum source term
     */
    void
    apply_free_atmosphere_damping(int lev, amrex::MultiFab& src_term) const;

    [[nodiscard]] amrex::Real forcing_height() const
    {
        return m_forcing_height;
    }

private:
    //! Read the free-atmosphere damping inputs
    void read_free_atmosphere_inputs(const amrex::ParmParse& pp_abl);

    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

    //! Damp the free atmosphere toward the geostrophic wind
    bool m_fa_damping{false};

    //! Set the free-atmosphere height from the capping inversion height
    bool m_fa_detect_height{false};

    //! Height above which the free atmosphere is damped
    amrex::Real m_fa_height{0.0_rt};

    //! Relaxation time scale of the free-atmosphere damping (s)
    amrex::Real m_fa_tau{100.0_rt};

    //! Time at which the free-atmosphere damping starts
    amrex::Real m_fa_start_time{0.0_rt};

    //! Time at which the free-atmosphere damping ends
    amrex::Real m_fa_end_time{std::numeric_limits<amrex::Real>::max()};

    //! Coriolis parameter f = 2 Omega sin(latitude)
    amrex::Real m_fa_coriolis_factor{1.0_rt};

    //! Number of free-atmosphere updates since start-up
    int m_fa_num_updates{0};

    //! Averaging direction of the velocity profile
    int m_fa_axis{2};

    //! Heights of the planar-averaged velocity profile
    amrex::Gpu::DeviceVector<amrex::Real> m_fa_heights;

    //! Planar-averaged velocity profile (3 components per height)
    amrex::Gpu::DeviceVector<amrex::Real> m_fa_velocity;

    //! Activated when water is present in domain
    bool m_use_phase_ramp{false};

    //! Number of cells in band to prevent forcing near liquid
    int m_n_band{2};

    //! ABL forcing terms
    amrex::RealArray m_abl_forcing{0.0_rt};

    //! File name for velocity forcing time table
    std::string m_vel_timetable;

    //! Bool for writing forcing time table
    bool m_write_force_timetable{false};
    //! File name for forcing time table output
    std::string m_force_timetable;
    //! Output frequency for forces
    int m_force_outfreq{1};
    //! Output start time for force
    amrex::Real m_force_outstart{0.0_rt};

    //! Velocity forcing time table
    amrex::Vector<amrex::Real> m_time_table;

    //! Velocity forcing speed table
    amrex::Vector<amrex::Real> m_speed_table;

    //! Velocity forcing direction table
    amrex::Vector<amrex::Real> m_direction_table;

    //! Target velocity
    amrex::Vector<amrex::Real> m_target_vel{0.0_rt, 0.0_rt, 0.0_rt};

    //! Current mean vel
    amrex::RealArray m_mean_vel{0.0_rt};

    //! Height at which the velocities are forcing
    amrex::Real m_forcing_height;

    //! Height from water interface where force is off
    amrex::Real m_forcing_mphase0;
    //! Height from water interface over which force is ramped
    amrex::Real m_forcing_mphase1;

    //! Local storage of interface location
    amrex::Real m_water_level;

    //! VOF field, to avoid forcing on liquid above force-off height
    const Field* m_vof;
};

} // namespace kynema_sgf::pde::icns

#endif /* ABLFORCING_H */
