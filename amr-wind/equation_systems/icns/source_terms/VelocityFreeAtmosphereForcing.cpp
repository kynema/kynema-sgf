#include "amr-wind/equation_systems/icns/source_terms/VelocityFreeAtmosphereForcing.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

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
    const auto& velocity =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const amrex::Real sponge_start = m_meso_start;
    const amrex::Real meso_timescale = m_time.delta_t();
    const auto vsize = m_wind_heights_d.size();
    const auto* wind_heights_d = m_wind_heights_d.data();
    const auto* u_values_d = m_u_values_d.data();
    const auto* v_values_d = m_v_values_d.data();
    const auto* w_values_d = m_w_values_d.data();
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
        amrex::Real ref_windx = velocity(i, j, k, 0);
        amrex::Real ref_windy = velocity(i, j, k, 1);
        amrex::Real ref_windz = velocity(i, j, k, 2);
        if (zi > 0) {
            ref_windx = (vsize > 0) ? interp::linear(
                                          wind_heights_d,
                                          wind_heights_d + vsize, u_values_d, z)
                                    : velocity(i, j, k, 0);
            ref_windy = (vsize > 0) ? interp::linear(
                                          wind_heights_d,
                                          wind_heights_d + vsize, v_values_d, z)
                                    : velocity(i, j, k, 1);
            ref_windz = (vsize > 0) ? interp::linear(
                                          wind_heights_d,
                                          wind_heights_d + vsize, w_values_d, z)
                                    : velocity(i, j, k, 2);
        }
        src_term(i, j, k, 0) -=
            1.0_rt / meso_timescale * (velocity(i, j, k, 0) - ref_windx);
        src_term(i, j, k, 1) -=
            1.0_rt / meso_timescale * (velocity(i, j, k, 1) - ref_windy);
        src_term(i, j, k, 2) -=
            1.0_rt / meso_timescale * (velocity(i, j, k, 2) - ref_windz);
    });
}

} // namespace amr_wind::pde::icns
