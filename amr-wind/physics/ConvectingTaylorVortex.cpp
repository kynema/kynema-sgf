#include <numbers>
#include "amr-wind/physics/ConvectingTaylorVortex.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ctv {

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real UExact::operator()(
    const amrex::Real u0,
    const amrex::Real v0,
    const amrex::Real omega,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real t) const
{
    return u0 - (std::cos(std::numbers::pi_v<amrex::Real> * (x - u0 * t)) *
                 std::sin(std::numbers::pi_v<amrex::Real> * (y - v0 * t)) *
                 std::exp(-2.0_rt * omega * t));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real VExact::operator()(
    const amrex::Real u0,
    const amrex::Real v0,
    const amrex::Real omega,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real t) const
{
    return v0 + (std::sin(std::numbers::pi_v<amrex::Real> * (x - u0 * t)) *
                 std::cos(std::numbers::pi_v<amrex::Real> * (y - v0 * t)) *
                 std::exp(-2.0_rt * omega * t));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real WExact::operator()(
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    return 0.0_rt;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real GpxExact::operator()(
    const amrex::Real u0,
    const amrex::Real /*unused*/,
    const amrex::Real omega,
    const amrex::Real x,
    const amrex::Real /*unused*/,
    const amrex::Real t) const
{
    return 0.5_rt * std::numbers::pi_v<amrex::Real> *
           std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * (x - u0 * t)) *
           std::exp(-4.0_rt * omega * t);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real GpyExact::operator()(
    const amrex::Real /*unused*/,
    const amrex::Real v0,
    const amrex::Real omega,
    const amrex::Real /*unused*/,
    const amrex::Real y,
    const amrex::Real t) const
{
    return 0.5_rt * std::numbers::pi_v<amrex::Real> *
           std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * (y - v0 * t)) *
           std::exp(-4.0_rt * omega * t);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real GpzExact::operator()(
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/,
    const amrex::Real /*unused*/) const
{
    return 0.0_rt;
}

} // namespace

ConvectingTaylorVortex::ConvectingTaylorVortex(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_gradp(sim.repo().get_field("gp"))
    , m_density(sim.repo().get_field("density"))
    , m_mesh_mapping(sim.has_mesh_mapping())
{
    {
        amrex::ParmParse pp("CTV");
        pp.query("density", m_rho);
        pp.query("u0", m_u0);
        pp.query("v0", m_v0);
        pp.query("activate_pressure", m_activate_pressure);
        pp.query("error_log_file", m_output_fname);
    }
    {
        amrex::Real nu;
        amrex::ParmParse pp("transport");
        pp.query("viscosity", nu);
        m_omega = std::numbers::pi_v<amrex::Real> *
                  std::numbers::pi_v<amrex::Real> * nu;
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u"
          << std::setw(m_w) << "L2_v" << std::setw(m_w) << "L2_w"
          << std::setw(m_w) << "L2_gpx" << std::setw(m_w) << "L2_gpy"
          << std::setw(m_w) << "L2_gpz" << '\n';
        f.close();
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void ConvectingTaylorVortex::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();

    const auto u0 = m_u0;
    const auto v0 = m_v0;
    const auto omega = m_omega;
    const bool activate_pressure = m_activate_pressure;
    const auto mesh_mapping = m_mesh_mapping;

    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& pressure = m_repo.get_field("p")(level);
    auto& gradp = m_repo.get_field("gp")(level);
    Field const* nu_coord_cc =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    Field const* nu_coord_nd =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_nd")) : nullptr;

    density.setVal(m_rho);

    UExact u_exact;
    VExact v_exact;
    WExact w_exact;
    GpxExact gpx_exact;
    GpyExact gpy_exact;
    GpzExact gpz_exact;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.tilebox();

        auto vel = velocity.array(mfi);
        auto gp = gradp.array(mfi);
        amrex::Array4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(level).const_array(mfi))
                         : amrex::Array4<amrex::Real const>();

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real x = mesh_mapping
                                ? (nu_cc(i, j, k, 0))
                                : (prob_lo[0] + ((i + 0.5_rt) * dx[0]));
            amrex::Real y = mesh_mapping
                                ? (nu_cc(i, j, k, 1))
                                : (prob_lo[1] + ((j + 0.5_rt) * dx[1]));

            vel(i, j, k, 0) = u_exact(u0, v0, omega, x, y, 0.0_rt);
            vel(i, j, k, 1) = v_exact(u0, v0, omega, x, y, 0.0_rt);
            vel(i, j, k, 2) = w_exact(u0, v0, omega, x, y, 0.0_rt);

            if (activate_pressure) {
                gp(i, j, k, 0) = gpx_exact(u0, v0, omega, x, y, 0.0_rt);
                gp(i, j, k, 1) = gpy_exact(u0, v0, omega, x, y, 0.0_rt);
                gp(i, j, k, 2) = gpz_exact(u0, v0, omega, x, y, 0.0_rt);
            }
        });

        if (activate_pressure) {
            const auto& nbx = mfi.nodaltilebox();
            auto pres = pressure.array(mfi);
            amrex::Array4<amrex::Real const> nu_nd =
                mesh_mapping ? ((*nu_coord_nd)(level).const_array(mfi))
                             : amrex::Array4<amrex::Real const>();

            amrex::ParallelFor(nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real x = mesh_mapping ? (nu_nd(i, j, k, 0))
                                             : (prob_lo[0] + (i * dx[0]));
                amrex::Real y = mesh_mapping ? (nu_nd(i, j, k, 1))
                                             : (prob_lo[1] + (j * dx[1]));

                pres(i, j, k, 0) =
                    -0.25_rt *
                    (std::cos(2.0_rt * std::numbers::pi_v<amrex::Real> * x) +
                     std::cos(2.0_rt * std::numbers::pi_v<amrex::Real> * y));
            });
        }
    }
}

template <typename T>
amrex::Real ConvectingTaylorVortex::compute_error(const Field& field)
{
    amrex::Real error = 0.0_rt;
    amrex::Real vol = 0.0_rt;
    const amrex::Real time = m_time.new_time();
    const auto u0 = m_u0;
    const auto v0 = m_v0;
    const auto omega = m_omega;
    T f_exact;
    const auto comp = f_exact.m_comp;
    const auto mesh_mapping = m_mesh_mapping;

    Field const* nu_coord_cc =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    Field const* mesh_fac_cc =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        if (m_sim.has_overset()) {
            const auto& iblank_arrs =
                m_repo.get_int_field("iblank_cell")(lev).const_arrays();
            const auto& imask_arrs = level_mask.arrays();
            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    if (std::abs(iblank_arrs[nbx](i, j, k)) < 1) {
                        imask_arrs[nbx](i, j, k) = 0;
                    }
                });
            amrex::Gpu::streamSynchronize();
        }

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& prob_lo = m_mesh.Geom(lev).ProbLoArray();

        const auto& fld = field(lev);
        auto const& fld_arr = fld.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();
        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();
        amrex::MultiArray4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        amrex::GpuTuple<amrex::Real, amrex::Real> vals = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real, amrex::Real>{}, fld, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real, amrex::Real> {
                auto const& fld_bx = fld_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                amrex::Real x = mesh_mapping
                                    ? (nu_cc[box_no](i, j, k, 0))
                                    : (prob_lo[0] + ((i + 0.5_rt) * dx[0]));
                amrex::Real y = mesh_mapping
                                    ? (nu_cc[box_no](i, j, k, 1))
                                    : (prob_lo[1] + ((j + 0.5_rt) * dx[1]));
                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0_rt;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0_rt;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0_rt;

                const amrex::Real u = fld_bx(i, j, k, comp);
                const amrex::Real u_exact = f_exact(u0, v0, omega, x, y, time);
                const amrex::Real cell_vol =
                    dx[0] * fac_x * dx[1] * fac_y * dx[2] * fac_z;

                return {
                    cell_vol * mask_bx(i, j, k) * (u - u_exact) * (u - u_exact),
                    cell_vol * mask_bx(i, j, k)};
            });
        error += amrex::get<0>(vals);
        vol += amrex::get<1>(vals);
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);
    amrex::ParallelDescriptor::ReduceRealSum(vol);

    return std::sqrt(error / vol);
}

void ConvectingTaylorVortex::output_error()
{
    // TODO: gradp analytical solution has not been adjusted for mesh mapping
    const amrex::Real u_err = compute_error<UExact>(m_velocity);
    const amrex::Real v_err = compute_error<VExact>(m_velocity);
    const amrex::Real w_err = compute_error<WExact>(m_velocity);
    const amrex::Real gpx_err = compute_error<GpxExact>(m_gradp);
    const amrex::Real gpy_err = compute_error<GpyExact>(m_gradp);
    const amrex::Real gpz_err = compute_error<GpzExact>(m_gradp);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time()
          << std::setw(m_w) << u_err << std::setw(m_w) << v_err
          << std::setw(m_w) << w_err << std::setw(m_w) << gpx_err
          << std::setw(m_w) << gpy_err << std::setw(m_w) << gpz_err << '\n';
        f.close();
    }
}

void ConvectingTaylorVortex::post_init_actions() { output_error(); }

void ConvectingTaylorVortex::post_advance_work() { output_error(); }

} // namespace amr_wind::ctv
