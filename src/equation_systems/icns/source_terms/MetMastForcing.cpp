#include "amr-wind/equation_systems/icns/source_terms/MetMastForcing.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

MetMastForcing::MetMastForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_sim(sim)
{
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("metmast_1dprofile_file", m_1d_metmast);
    if (!m_1d_metmast.empty()) {
        std::ifstream ransfile(m_1d_metmast, std::ios::in);
        if (!ransfile.good()) {
            amrex::Abort("Cannot find Met Mast profile file " + m_1d_metmast);
        }
        //! x y z u v w T
        amrex::Real value1, value2, value3, value4, value5, value6, value7;
        while (ransfile >> value1 >> value2 >> value3 >> value4 >> value5 >>
               value6 >> value7) {
            m_metmast_x.push_back(value1);
            m_metmast_y.push_back(value2);
            m_metmast_z.push_back(value3);
            m_u_values.push_back(value4);
            m_v_values.push_back(value5);
            m_w_values.push_back(value6);
        }
    } else {
        amrex::Abort("Cannot find 1-D Met Mast profile file " + m_1d_metmast);
    }
    pp_abl.query("meso_timescale", m_meso_timescale);
    pp_abl.query("metmast_horizontal_radius", m_long_radius);
    pp_abl.query("metmast_vertical_radius", m_vertical_radius);
    pp_abl.query("metmast_damping_radius", m_damping_radius);
    int num_wind_values = static_cast<int>(m_metmast_x.size());
    m_metmast_x_d.resize(num_wind_values);
    m_metmast_y_d.resize(num_wind_values);
    m_metmast_z_d.resize(num_wind_values);
    m_u_values_d.resize(num_wind_values);
    m_v_values_d.resize(num_wind_values);
    m_w_values_d.resize(num_wind_values);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_x.begin(), m_metmast_x.end(),
        m_metmast_x_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_y.begin(), m_metmast_y.end(),
        m_metmast_y_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_metmast_z.begin(), m_metmast_z.end(),
        m_metmast_z_d.begin());
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

MetMastForcing::~MetMastForcing() = default;

void MetMastForcing::operator()(
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
    const amrex::Real meso_timescale = m_meso_timescale;
    const amrex::Real long_radius = m_long_radius;
    const amrex::Real vertical_radius = m_vertical_radius;
    const amrex::Real damping_radius = m_damping_radius;
    const auto vsize = m_metmast_x_d.size();
    const auto* metmast_x_d = m_metmast_x_d.data();
    const auto* metmast_y_d = m_metmast_y_d.data();
    const auto* metmast_z_d = m_metmast_z_d.data();
    const auto* u_values_d = m_u_values_d.data();
    const auto* v_values_d = m_v_values_d.data();
    const auto* w_values_d = m_w_values_d.data();
    const bool has_terrain = this->m_sim.repo().field_exists("terrain_height");
    if (has_terrain) {
        auto* const m_terrain_height =
            &this->m_sim.repo().get_field("terrain_height");
        const auto& terrain_height = (*m_terrain_height)(lev).const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = amrex::max<amrex::Real>(
                prob_lo[2] + (k + 0.5_rt) * dx[2] - terrain_height(i, j, k),
                0.5_rt * dx[2]);
            //! Testing for one point for now
            int ii = 0;
            amrex::Real ri2 = (x - metmast_x_d[ii]) * (x - metmast_x_d[ii]) /
                              (long_radius * long_radius);
            ri2 += (y - metmast_y_d[ii]) * (y - metmast_y_d[ii]) /
                   (long_radius * long_radius);
            ri2 += (z - metmast_z_d[ii]) * (z - metmast_z_d[ii]) /
                   (vertical_radius * vertical_radius);
            const amrex::Real weight_fn =
                (ri2 <= damping_radius) ? std::exp(-0.25_rt * ri2) : 0;
            amrex::Real ref_windx = u_values_d[ii];
            amrex::Real ref_windy = v_values_d[ii];
            amrex::Real ref_windz = w_values_d[ii];
            src_term(i, j, k, 0) -=
                weight_fn / meso_timescale * (velocity(i, j, k, 0) - ref_windx);
            src_term(i, j, k, 1) -=
                weight_fn / meso_timescale * (velocity(i, j, k, 1) - ref_windy);
            src_term(i, j, k, 2) -=
                weight_fn / meso_timescale * (velocity(i, j, k, 2) - ref_windz);
        });
    } else {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5_rt) * dx[2];
            //! Testing for one point for now
            int ii = 0;
            amrex::Real ri2 = (x - metmast_x_d[ii]) * (x - metmast_x_d[ii]) /
                              (long_radius * long_radius);
            ri2 += (y - metmast_y_d[ii]) * (y - metmast_y_d[ii]) /
                   (long_radius * long_radius);
            ri2 += (z - metmast_z_d[ii]) * (z - metmast_z_d[ii]) /
                   (vertical_radius * vertical_radius);
            const amrex::Real weight_fn =
                (ri2 <= 100) ? std::exp(-0.25_rt * ri2) : 0;
            amrex::Real ref_windx = u_values_d[ii];
            amrex::Real ref_windy = v_values_d[ii];
            amrex::Real ref_windz = w_values_d[ii];
            src_term(i, j, k, 0) -=
                weight_fn / meso_timescale * (velocity(i, j, k, 0) - ref_windx);
            src_term(i, j, k, 1) -=
                weight_fn / meso_timescale * (velocity(i, j, k, 1) - ref_windy);
            src_term(i, j, k, 2) -=
                weight_fn / meso_timescale * (velocity(i, j, k, 2) - ref_windz);
        });
    }
}

} // namespace amr_wind::pde::icns
