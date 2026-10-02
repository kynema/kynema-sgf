#include "src/equation_systems/icns/source_terms/ABLForcing.H"
#include "src/CFDSim.H"
#include "src/wind_energy/ABL.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "src/equation_systems/vof/volume_fractions.H"
#include "src/utilities/trig_ops.H"

#include <algorithm>

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

ABLForcing::ABLForcing(const CFDSim& sim)
    : m_time(sim.time()), m_mesh(sim.mesh())
{
    const auto& abl = sim.physics_manager().get<kynema_sgf::ABL>();
    abl.register_forcing_term(this);
    abl.abl_statistics().register_forcing_term(this);

    amrex::ParmParse pp_abl(identifier());
    // TODO: Allow forcing at multiple heights
    pp_abl.get("abl_forcing_height", m_forcing_height);
    amrex::ParmParse pp_incflo("incflo");

    pp_abl.query("velocity_timetable", m_vel_timetable);
    if (!m_vel_timetable.empty()) {
        std::ifstream ifh(m_vel_timetable, std::ios::in);
        if (!ifh.good()) {
            amrex::Abort(
                "Cannot find ABLForcing velocity_timetable file: " +
                m_vel_timetable);
        }
        amrex::Real data_time;
        amrex::Real data_speed;
        amrex::Real data_deg;
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        while (ifh >> data_time) {
            ifh >> data_speed >> data_deg;
            amrex::Real data_rad = utils::radians(data_deg);
            m_time_table.push_back(data_time);
            m_speed_table.push_back(data_speed);
            m_direction_table.push_back(data_rad);
        }
    } else {
        pp_incflo.getarr("velocity", m_target_vel);
    }

    pp_abl.query("free_atmosphere_damping", m_fa_damping);
    if (m_fa_damping) {
        read_free_atmosphere_inputs(pp_abl);
    }

    m_write_force_timetable = pp_abl.contains("forcing_timetable_output_file");
    if (m_write_force_timetable) {
        pp_abl.get("forcing_timetable_output_file", m_force_timetable);
        pp_abl.query("forcing_timetable_frequency", m_force_outfreq);
        pp_abl.query("forcing_timetable_start_time", m_force_outstart);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream outfile;
            outfile.open(m_force_timetable, std::ios::out);
            outfile << "time\tfx\tfy\tfz";
            if (m_fa_damping) {
                outfile << "\tugx\tugy\tfa_height";
            }
            outfile << '\n';
        }
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_mean_vel[i] = m_target_vel[i];
    }

    // Set up relaxation toward 0 forcing near the air-water interface
    if (sim.repo().field_exists("vof")) {
        // If vof exists, get multiphase physics
        const auto& mphase =
            sim.physics_manager().get<kynema_sgf::MultiPhase>();
        // Retrieve interface position
        m_water_level = mphase.water_level();
        // Confirm that water level will be used
        m_use_phase_ramp = true;
        // Parse for thickness of ramping function
        pp_abl.get("abl_forcing_off_height", m_forcing_mphase0);
        pp_abl.get("abl_forcing_ramp_height", m_forcing_mphase1);
        // Store reference to vof field
        m_vof = &sim.repo().get_field("vof");
        // Parse for number of cells in band
        pp_abl.query("abl_forcing_band", m_n_band);
    } else {
        // Point to something, will not be used
        m_vof = &sim.repo().get_field("velocity");
    }
}

ABLForcing::~ABLForcing() = default;

void ABLForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    const amrex::Real dudt = m_abl_forcing[0];
    const amrex::Real dvdt = m_abl_forcing[1];

    const bool ph_ramp = m_use_phase_ramp;
    const int n_band = m_n_band;
    const amrex::Real wlev = m_water_level;
    const amrex::Real wrht0 = m_forcing_mphase0;
    const amrex::Real wrht1 = m_forcing_mphase1;
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    auto const& src_arrs = src_term.arrays();
    auto const& vof_arrs = (*m_vof)(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            if (n >= 2) {
                return;
            }
            amrex::Real fac = 1.0_rt;
            if (ph_ramp) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                if (z - wlev < wrht0 + wrht1) {
                    if (z - wlev < wrht0) {
                        fac = 0.0_rt;
                    } else {
                        fac = 0.5_rt -
                              (0.5_rt * std::cos(
                                            std::numbers::pi_v<amrex::Real> *
                                            (z - wlev - wrht0) / wrht1));
                    }
                }
                if (multiphase::interface_band(
                        i, j, k, vof_arrs[nbx], n_band) ||
                    vof_arrs[nbx](i, j, k) >
                        1.0_rt - (std::numeric_limits<amrex::Real>::epsilon() *
                                  1.0e4_rt)) {
                    fac = 0.0_rt;
                }
            }
            const amrex::Real forcing = (n == 0) ? dudt : dvdt;
            src_arrs[nbx](i, j, k, n) += fac * forcing;
        });

    if (m_fa_damping) {
        apply_free_atmosphere_damping(lev, src_term);
    }
}

void ABLForcing::read_free_atmosphere_inputs(const amrex::ParmParse& pp_abl)
{
    pp_abl.query("detect_free_atmosphere_height", m_fa_detect_height);
    if (!m_fa_detect_height) {
        pp_abl.get("free_atmosphere_height", m_fa_height);
    }
    pp_abl.query("free_atmosphere_damping_time_scale", m_fa_tau);
    pp_abl.query("free_atmosphere_damping_start_time", m_fa_start_time);
    pp_abl.query("free_atmosphere_damping_end_time", m_fa_end_time);
    if (m_fa_tau <= 0.0_rt) {
        amrex::Abort(
            "ABLForcing: free_atmosphere_damping_time_scale must be positive");
    }

    // The geostrophic wind only balances the forcing when the Coriolis force
    // acts on the flow
    amrex::Vector<std::string> icns_sources;
    amrex::ParmParse("ICNS").queryarr("source_terms", icns_sources);
    if (std::ranges::find(icns_sources, "CoriolisForcing") ==
        icns_sources.end()) {
        amrex::Abort(
            "ABLForcing: free_atmosphere_damping requires CoriolisForcing in "
            "ICNS.source_terms");
    }

    // Same Coriolis parameter as GeostrophicForcing
    amrex::ParmParse ppc("CoriolisForcing");
    amrex::Real rot_time_period = 86164.091_rt;
    ppc.query("rotational_time_period", rot_time_period);
    amrex::Real latitude = 0.0_rt;
    ppc.get("latitude", latitude);
    const amrex::Real sinphi = std::sin(utils::radians(latitude));
    if (std::abs(sinphi) < 1.0e-6_rt) {
        amrex::Abort(
            "ABLForcing: free_atmosphere_damping needs a nonzero Coriolis "
            "parameter (latitude away from the equator)");
    }
    m_fa_coriolis_factor = 2.0_rt * utils::two_pi() / rot_time_period * sinphi;
}

void ABLForcing::update_free_atmosphere(
    const VelPlaneAveraging& vpa, const amrex::Real zi)
{
    if (!m_fa_damping) {
        return;
    }

    ++m_fa_num_updates;
    m_fa_axis = vpa.axis();
    if (m_fa_detect_height) {
        m_fa_height = m_mesh.Geom(0).ProbLo(m_fa_axis) + zi;
    }

    const auto& heights = vpa.line_centroids();
    const auto& velocity = vpa.line_average();
    m_fa_heights.resize(heights.size());
    m_fa_velocity.resize(velocity.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, heights.begin(), heights.end(),
        m_fa_heights.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, velocity.begin(), velocity.end(),
        m_fa_velocity.begin());
}

void ABLForcing::apply_free_atmosphere_damping(
    const int lev, amrex::MultiFab& src_term) const
{
    // The forcing at a timestep corrects the drift of the previous timestep,
    // so the first forcing after start-up (zero on a cold start) does not
    // carry the geostrophic wind yet. The initial iterations run before any
    // update.
    if (m_fa_num_updates < 2) {
        return;
    }
    const amrex::Real time =
        0.5_rt * (m_time.current_time() + m_time.new_time());
    if ((time < m_fa_start_time) || (time > m_fa_end_time)) {
        return;
    }

    const auto ug = geostrophic_velocity();
    const amrex::Real ugx = ug[0];
    const amrex::Real ugy = ug[1];
    const amrex::Real inv_tau = 1.0_rt / m_fa_tau;
    const amrex::Real fa_height = m_fa_height;
    const int idir = m_fa_axis;
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    const amrex::Real* heights = m_fa_heights.data();
    const amrex::Real* heights_end = m_fa_heights.end();
    const amrex::Real* vel = m_fa_velocity.data();

    auto const& src_arrs = src_term.arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const amrex::Real z =
                problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);
            if (z < fa_height) {
                return;
            }
            const amrex::Real umean = kynema_sgf::interp::linear(
                heights, heights_end, vel, z, AMREX_SPACEDIM, 0);
            const amrex::Real vmean = kynema_sgf::interp::linear(
                heights, heights_end, vel, z, AMREX_SPACEDIM, 1);
            src_arrs[nbx](i, j, k, 0) += inv_tau * (ugx - umean);
            src_arrs[nbx](i, j, k, 1) += inv_tau * (ugy - vmean);
        });
}

} // namespace kynema_sgf::pde::icns
