#include "src/equation_systems/icns/source_terms/VelocityFreeAtmosphereForcing.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

VelocityFreeAtmosphereForcing::VelocityFreeAtmosphereForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_sim(sim)
{
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("rans_1dprofile_file", m_1d_rans_filename);
    if (!m_1d_rans_filename.empty()) {
        std::ifstream ransfile(m_1d_rans_filename, std::ios::in);
        if (!ransfile.good()) {
            amrex::Abort(
                "Cannot find 1-D RANS profile file " + m_1d_rans_filename);
        }
        amrex::Real value1, value2, value3, value4, value5;
        while (ransfile >> value1 >> value2 >> value3 >> value4 >> value5) {
            m_wind_heights.push_back(value1);
            m_u_values.push_back(value2);
            m_v_values.push_back(value3);
            m_w_values.push_back(value4);
        }
    } else {
        amrex::Abort("Cannot find 1-D RANS profile file " + m_1d_rans_filename);
    }
    pp_abl.query("meso_sponge_start", m_meso_start);
    pp_abl.query("velocity_sponge_start", m_meso_start);
    pp_abl.query("meso_timescale", m_meso_timescale);
    int num_wind_values = static_cast<int>(m_wind_heights.size());
    m_wind_heights_d.resize(num_wind_values);
    m_u_values_d.resize(num_wind_values);
    m_v_values_d.resize(num_wind_values);
    m_w_values_d.resize(num_wind_values);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_wind_heights.begin(), m_wind_heights.end(),
        m_wind_heights_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_u_values.begin(), m_u_values.end(),
        m_u_values_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_v_values.begin(), m_v_values.end(),
        m_v_values_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_w_values.begin(), m_w_values.end(),
        m_w_values_d.begin());
}

VelocityFreeAtmosphereForcing::~VelocityFreeAtmosphereForcing() = default;

void VelocityFreeAtmosphereForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();

    const amrex::Real sponge_start = m_meso_start;
    const amrex::Real meso_timescale = m_time.delta_t();
    const auto vsize = m_wind_heights_d.size();
    const auto* wind_heights_d = m_wind_heights_d.data();
    const auto* u_values_d = m_u_values_d.data();
    const auto* v_values_d = m_v_values_d.data();
    const auto* w_values_d = m_w_values_d.data();

    const bool has_terrain = this->m_sim.repo().field_exists("terrain_height");
    auto const& terrain_arrs =
        has_terrain
            ? this->m_sim.repo().get_field("terrain_height")(lev).const_arrays()
            : amrex::MultiArray4<amrex::Real const>();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            const amrex::Real cell_terrain_height =
                (has_terrain) ? terrain_arrs[nbx](i, j, k) : 0.0_rt;
            const amrex::Real z = amrex::max<amrex::Real>(
                prob_lo[2] + ((k + 0.5_rt) * dx[2]) - cell_terrain_height,
                0.5_rt * dx[2]);
            const amrex::Real zi = amrex::max<amrex::Real>(
                (z - sponge_start) / (prob_hi[2] - sponge_start), 0.0_rt);
            amrex::Real ref_wind = vel_arrs[nbx](i, j, k, n);
            if (zi > 0 && vsize > 0) {
                const amrex::Real* wind_ht_end = wind_heights_d + vsize;
                const amrex::Real* vals = (n == 0)   ? u_values_d
                                          : (n == 1) ? v_values_d
                                                     : w_values_d;
                ref_wind = kynema_sgf::interp::linear(
                    wind_heights_d, wind_ht_end, vals, z);
            }
            src_arrs[nbx](i, j, k, n) -= 1.0_rt / meso_timescale *
                                         (vel_arrs[nbx](i, j, k, n) - ref_wind);
        });
}

} // namespace kynema_sgf::pde::icns
