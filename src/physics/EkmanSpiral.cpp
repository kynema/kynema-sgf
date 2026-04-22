#include "amr-wind/physics/EkmanSpiral.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

namespace {

struct UExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real /*v0*/, amrex::Real /*a*/, amrex::Real /*z*/) const;
    const int m_comp{0};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real UExact::operator()(
    const amrex::Real v0, const amrex::Real a, const amrex::Real z) const
{
    return v0 * (1.0_rt - std::exp(-a * z) * std::cos(-a * z));
}

struct VExact
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real /*v0*/, amrex::Real /*a*/, amrex::Real /*z*/) const;
    const int m_comp{1};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real VExact::operator()(
    const amrex::Real v0, const amrex::Real a, const amrex::Real z) const
{
    return -v0 * std::exp(-a * z) * std::sin(-a * z);
}

} // namespace

EkmanSpiral::EkmanSpiral(const CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{

    {
        amrex::ParmParse pp("incflo");
        pp.query("density", m_rho);
    }

    amrex::Real coriolis_factor;
    {
        amrex::ParmParse pp("CoriolisForcing");
        amrex::Real rot_time_period;
        pp.get("rotational_time_period", rot_time_period);
        coriolis_factor = 2.0_rt * utils::two_pi() / rot_time_period;

        amrex::Real latitude;
        pp.get("latitude", latitude);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs(latitude - 90.0_rt) <
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt,
            "Ekman Spiral only works with geostrophic forcing which has to be "
            "at latitude 90 degrees");
    }

    {
        amrex::ParmParse pp("GeostrophicForcing");
        amrex::Vector<amrex::Real> gwind{15.0_rt, 0.0_rt, 0.0_rt};
        pp.getarr("geostrophic_wind", gwind);
        m_vel = gwind[0];

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs(gwind[1]) <
                    std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt &&
                std::abs(gwind[2]) <
                    std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt,
            "Ekman Spiral only works for forcing in x-dir for now");
    }

    amrex::Real Az;
    {
        amrex::ParmParse pp("transport");
        pp.get("viscosity", Az);
    }

    m_DE = std::sqrt(2.0_rt * Az / coriolis_factor);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u"
          << std::setw(m_w) << "L2_v" << '\n';
        f.close();
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void EkmanSpiral::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    density.setVal(m_rho);

    amrex::Real a = 1.0_rt / m_DE;
    amrex::Real v0 = m_vel;

    UExact u_exact;
    VExact v_exact;

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& vel_arrs = velocity.arrays();

    amrex::ParallelFor(
        velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            vel_arrs[nbx](i, j, k, 0) = u_exact(v0, a, z);
            vel_arrs[nbx](i, j, k, 1) = v_exact(v0, a, z);
            vel_arrs[nbx](i, j, k, 2) = 0.0_rt;
        });
    amrex::Gpu::streamSynchronize();
}

template <typename T>
amrex::Real EkmanSpiral::compute_error(const Field& field)
{

    amrex::Real error = 0.0_rt;
    const auto v0 = m_vel;
    const auto a = 1.0_rt / m_DE;

    T f_exact;
    const auto comp = f_exact.m_comp;

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

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& problo = m_mesh.Geom(lev).ProbLoArray();
        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        const auto& fld = field(lev);
        error += amrex::ReduceSum(
            fld, level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& fld_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real err_fab = 0.0_rt;

                amrex::Loop(bx, [=, &err_fab](int i, int j, int k) {
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                    const amrex::Real u = fld_arr(i, j, k, comp);
                    const amrex::Real u_exact = f_exact(v0, a, z);
                    err_fab += cell_vol * mask_arr(i, j, k) * (u - u_exact) *
                               (u - u_exact);
                });
                return err_fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void EkmanSpiral::output_error()
{
    const amrex::Real u_err = compute_error<UExact>(m_velocity);
    const amrex::Real v_err = compute_error<VExact>(m_velocity);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time()
          << std::setw(m_w) << u_err << std::setw(m_w) << v_err << '\n';
        f.close();
    }
}

void EkmanSpiral::post_init_actions() { output_error(); }

void EkmanSpiral::post_advance_work() { output_error(); }

} // namespace amr_wind
