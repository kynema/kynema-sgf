#include "src/equation_systems/icns/source_terms/MetMastForcing.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

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
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();

    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();

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
    auto const& terrain_arrs =
        has_terrain
            ? this->m_sim.repo().get_field("terrain_height")(lev).const_arrays()
            : amrex::MultiArray4<amrex::Real const>();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            const amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = has_terrain
                                      ? amrex::max<amrex::Real>(
                                            prob_lo[2] + (k + 0.5_rt) * dx[2] -
                                                terrain_arrs[nbx](i, j, k),
                                            0.5_rt * dx[2])
                                      : prob_lo[2] + (k + 0.5_rt) * dx[2];
            int ii = 0;
            amrex::Real ri2 = (x - metmast_x_d[ii]) * (x - metmast_x_d[ii]) /
                              (long_radius * long_radius);
            ri2 += (y - metmast_y_d[ii]) * (y - metmast_y_d[ii]) /
                   (long_radius * long_radius);
            ri2 += (z - metmast_z_d[ii]) * (z - metmast_z_d[ii]) /
                   (vertical_radius * vertical_radius);
            const amrex::Real weight_fn =
                (ri2 <= (has_terrain ? damping_radius : 100.0_rt))
                    ? std::exp(-0.25_rt * ri2)
                    : 0;
            const amrex::Real* vals = (n == 0)   ? u_values_d
                                      : (n == 1) ? v_values_d
                                                 : w_values_d;
            const amrex::Real ref_wind = vals[ii];
            src_arrs[nbx](i, j, k, n) -= weight_fn / meso_timescale *
                                         (vel_arrs[nbx](i, j, k, n) - ref_wind);
        });
}

} // namespace kynema_sgf::pde::icns
