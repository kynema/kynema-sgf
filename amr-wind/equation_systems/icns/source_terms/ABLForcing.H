#ifndef ABLFORCING_H
#define ABLFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Forcing term to drive the ABL profile to a desired state
 *  \ingroup icns_src we_abl
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
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

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
            const amrex::Real new_spd = ::amr_wind::interp::linear(
                m_time_table, m_speed_table, new_time);
            const amrex::Real new_dir = ::amr_wind::interp::linear_angle(
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
                    << 0.0_rt << '\n';
        }
    }

    [[nodiscard]] amrex::RealArray abl_forcing() const { return m_abl_forcing; }

    [[nodiscard]] amrex::Real forcing_height() const
    {
        return m_forcing_height;
    }

private:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

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

} // namespace amr_wind::pde::icns

#endif /* ABLFORCING_H */
