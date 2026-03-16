#include "amr-wind/equation_systems/temperature/source_terms/TemperatureFreeAtmosphereForcing.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

TemperatureFreeAtmosphereForcing::TemperatureFreeAtmosphereForcing(
    const CFDSim& sim)
    : m_mesh(sim.mesh())
    , m_temperature(sim.repo().get_field("temperature"))
    , m_sim(sim)
{
    amrex::ParmParse pp_abl("ABL");
    //! Temperature variation as a function of height
    pp_abl.query("meso_sponge_start", m_meso_start);
    pp_abl.query("meso_timescale", m_meso_timescale);
    pp_abl.getarr("temperature_heights", m_theta_heights);
    pp_abl.getarr("temperature_values", m_theta_values);
    pp_abl.query("horizontal_sponge_temp", m_horizontal_sponge);
    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    const int num_theta_values = static_cast<int>(m_theta_heights.size());
    m_theta_heights_d.resize(num_theta_values);
    m_theta_values_d.resize(num_theta_values);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_heights.begin(),
        m_theta_heights.end(), m_theta_heights_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
        m_theta_values_d.begin());
    amrex::ParmParse pp("DragForcing");
    pp.query("sponge_strength", m_sponge_strength);
    pp.query("sponge_density", m_sponge_density);
    pp.query("sponge_distance_west", m_sponge_distance_west);
    pp.query("sponge_distance_east", m_sponge_distance_east);
    pp.query("sponge_distance_south", m_sponge_distance_south);
    pp.query("sponge_distance_north", m_sponge_distance_north);
    pp.query("sponge_west", m_sponge_west);
    pp.query("sponge_east", m_sponge_east);
    pp.query("sponge_south", m_sponge_south);
    pp.query("sponge_north", m_sponge_north);
}

TemperatureFreeAtmosphereForcing::~TemperatureFreeAtmosphereForcing() = default;

void TemperatureFreeAtmosphereForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const auto& temperature =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_array(
            mfi);
    const amrex::Real sponge_start = m_meso_start;
    const amrex::Real meso_timescale = m_meso_timescale;
    const auto vsize = m_theta_heights_d.size();
    const auto* theta_heights_d = m_theta_heights_d.data();
    const auto* theta_values_d = m_theta_values_d.data();
    const bool has_terrain = this->m_sim.repo().field_exists("terrain_height");
    auto* const m_terrain_height =
        (has_terrain) ? &this->m_sim.repo().get_field("terrain_height")
                      : nullptr;
    const auto& terrain_height = (has_terrain)
                                     ? (*m_terrain_height)(lev).const_array(mfi)
                                     : amrex::Array4<amrex::Real>();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real cell_terrain_height =
            (has_terrain) ? terrain_height(i, j, k) : 0.0_rt;
        const amrex::Real z = amrex::max<amrex::Real>(
            prob_lo[2] + ((k + 0.5_rt) * dx[2]) - cell_terrain_height,
            0.5_rt * dx[2]);
        const amrex::Real zi = amrex::max<amrex::Real>(
            (z - sponge_start) / (prob_hi[2] - sponge_start), 0.0_rt);
        amrex::Real ref_temp = temperature(i, j, k);
        if (zi > 0) {
            ref_temp = (vsize > 0)
                           ? interp::linear(
                                 theta_heights_d, theta_heights_d + vsize,
                                 theta_values_d, z)
                           : temperature(i, j, k);
        }
        src_term(i, j, k, 0) -=
            1.0_rt / meso_timescale * (temperature(i, j, k) - ref_temp);
    });
    if (m_horizontal_sponge) {
        const amrex::Real sponge_strength = m_sponge_strength;
        const amrex::Real sponge_density = m_sponge_density;
        const amrex::Real start_east = prob_hi[0] - m_sponge_distance_east;
        const amrex::Real start_west = prob_lo[0] - m_sponge_distance_west;
        const amrex::Real start_north = prob_hi[1] - m_sponge_distance_north;
        const amrex::Real start_south = prob_lo[1] - m_sponge_distance_south;
        const int sponge_east = m_sponge_east;
        const int sponge_west = m_sponge_west;
        const int sponge_south = m_sponge_south;
        const int sponge_north = m_sponge_north;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            amrex::Real xstart_damping = 0;
            amrex::Real ystart_damping = 0;
            amrex::Real xend_damping = 0;
            amrex::Real yend_damping = 0;
            amrex::Real xi_end = (x - start_east) / (prob_hi[0] - start_east);
            amrex::Real xi_start = (start_west - x) / (start_west - prob_lo[0]);
            xi_start = sponge_west * amrex::max<amrex::Real>(xi_start, 0.0_rt);
            xi_end = sponge_east * amrex::max<amrex::Real>(xi_end, 0.0_rt);
            xstart_damping =
                sponge_west * sponge_strength * xi_start * xi_start;
            xend_damping = sponge_east * sponge_strength * xi_end * xi_end;
            amrex::Real yi_end = (y - start_north) / (prob_hi[1] - start_north);
            amrex::Real yi_start =
                (start_south - y) / (start_south - prob_lo[1]);
            yi_start = sponge_south * amrex::max<amrex::Real>(yi_start, 0.0_rt);
            yi_end = sponge_north * amrex::max<amrex::Real>(yi_end, 0.0_rt);
            ystart_damping = sponge_strength * yi_start * yi_start;
            yend_damping = sponge_strength * yi_end * yi_end;
            const amrex::Real ref_temp =
                (vsize > 0) ? interp::linear(
                                  theta_heights_d, theta_heights_d + vsize,
                                  theta_values_d, z)
                            : temperature(i, j, k);
            src_term(i, j, k, 0) -=
                (xstart_damping + xend_damping + ystart_damping +
                 yend_damping) *
                (temperature(i, j, k) - sponge_density * ref_temp);
        });
    }
}

} // namespace amr_wind::pde::temperature
