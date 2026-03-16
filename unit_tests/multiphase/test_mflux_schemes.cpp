#include "aw_test_utils/MeshTest.H"
#include "hydro_godunov_ppm.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_scalar_increasing(amr_wind::Field& fld, int dir)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                farrs[nbx](i, j, k, 0) =
                    (amrex::Real)(dir == 0 ? i * i
                                           : (dir == 1 ? j * j : k * k));
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_scalar_slopechange(amr_wind::Field& fld, int dir, int center)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                farrs[nbx](i, j, k, 0) =
                    (amrex::Real)(dir == 0 ? std::abs(i - center)
                                           : (dir == 1 ? std::abs(j - center)
                                                       : std::abs(k - center)));
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_scalar_uniform(amr_wind::Field& fld, amrex::Real cst)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                farrs[nbx](i, j, k, 0) = cst;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void get_output_upwind(
    amr_wind::Field& fld,
    amr_wind::Field& mac_fld,
    amrex::Real dt,
    int ii,
    int jj,
    int kk,
    amrex::Real& Im,
    amrex::Real& Ip)
{
    const int nlevels = fld.repo().num_active_levels();
    amrex::Gpu::DeviceVector<amrex::Real> dout(2, 0.0_rt);
    auto* dout_ptr = dout.data();
    const auto* pbc = fld.bcrec_device().data();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = fld.repo().mesh().Geom(lev).CellSizeArray();
        amrex::Box const& domain = fld.repo().mesh().Geom(lev).Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        auto limiter = PPM::upwind();
        const auto& farrs = fld(lev).arrays();
        const auto& vel_mac_arrs = mac_fld(lev).const_arrays();

        amrex::ParallelFor(
            fld(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real im_tmp, ip_tmp;
                PPM::PredictStateOnXFace(
                    i, j, k, 0, dt, dx[0], im_tmp, ip_tmp, farrs[nbx],
                    vel_mac_arrs[nbx], pbc[0], dlo.x, dhi.x, limiter,
                    PPM::UPWIND);
                if (i == ii && j == jj && k == kk) {
                    dout_ptr[0] = im_tmp;
                    dout_ptr[1] = ip_tmp;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dout.begin(), dout.begin() + 1, &Im);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, dout.end() - 1, dout.end(), &Ip);
}

void get_output_minmod(
    amr_wind::Field& fld,
    amr_wind::Field& mac_fld,
    amrex::Real dt,
    int ii,
    int jj,
    int kk,
    int dir,
    amrex::Real& Im,
    amrex::Real& Ip)
{
    const int nlevels = fld.repo().num_active_levels();
    amrex::Gpu::DeviceVector<amrex::Real> dout(2, 0.0_rt);
    auto* dout_ptr = dout.data();
    const auto* pbc = fld.bcrec_device().data();

    auto limiter = PPM::minmod();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = fld.repo().mesh().Geom(lev).CellSizeArray();
        amrex::Box const& domain = fld.repo().mesh().Geom(lev).Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);
        const auto& farrs = fld(lev).arrays();
        const auto& vel_mac_arrs = mac_fld(lev).const_arrays();

        amrex::ParallelFor(
            fld(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real im_tmp = 0.0_rt;
                amrex::Real ip_tmp = 0.0_rt;
                if (dir == 0) {
                    PPM::PredictStateOnXFace(
                        i, j, k, 0, dt, dx[0], im_tmp, ip_tmp, farrs[nbx],
                        vel_mac_arrs[nbx], pbc[0], dlo.x, dhi.x, limiter,
                        PPM::MINMOD);
                } else if (dir == 1) {
                    PPM::PredictStateOnYFace(
                        i, j, k, 0, dt, dx[1], im_tmp, ip_tmp, farrs[nbx],
                        vel_mac_arrs[nbx], pbc[0], dlo.y, dhi.y, limiter,
                        PPM::MINMOD);
                } else if (dir == 2) {
                    PPM::PredictStateOnZFace(
                        i, j, k, 0, dt, dx[2], im_tmp, ip_tmp, farrs[nbx],
                        vel_mac_arrs[nbx], pbc[0], dlo.z, dhi.z, limiter,
                        PPM::MINMOD);
                }
                if (i == ii && j == jj && k == kk) {
                    dout_ptr[0] = im_tmp;
                    dout_ptr[1] = ip_tmp;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dout.begin(), dout.begin() + 1, &Im);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, dout.end() - 1, dout.end(), &Ip);
}

} // namespace

class MFluxSchemeTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", m_nx);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            pp.addarr("prob_lo", m_problo);
            pp.addarr("prob_hi", m_probhi);
        }
    }
    // Parameters
    const amrex::Vector<amrex::Real> m_problo{{0.0_rt, -1.0_rt, 0.0_rt}};
    const amrex::Vector<amrex::Real> m_probhi{{2.0_rt, +1.0_rt, 2.0_rt}};
    const int m_nx = 5;
    const amrex::Real m_tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
};

TEST_F(MFluxSchemeTest, upwind)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Initialize variables to be reused
    amrex::Real Im = 0.0_rt;
    amrex::Real Ip = 0.0_rt;

    // Initialize field variable
    auto& sc = repo.declare_field("scalar", 1, 2);

    // Parameters to be referenced during test
    amrex::Real adv_vel = 1.2_rt;
    amrex::Real dt = 0.25_rt;

    // Initialize mac velocity
    repo.declare_face_normal_field({"umac", "vmac", "wmac"}, 1, 1, 1);
    auto& umac = repo.get_field("umac");
    auto& vmac = repo.get_field("vmac");
    auto& wmac = repo.get_field("wmac");
    umac.setVal(adv_vel);
    vmac.setVal(adv_vel);
    wmac.setVal(-adv_vel);

    /* -- Constant field portion of test -- */
    // Set up field
    amrex::Real sc_cst = 1.5_rt;
    init_scalar_uniform(sc, sc_cst);
    // Compute interpolated quantities at each face
    int i = m_nx / 2;
    int j = m_nx / 2;
    int k = m_nx / 2;
    get_output_upwind(sc, umac, dt, i, j, k, Im, Ip);
    // Check values
    EXPECT_NEAR(sc_cst, Im, m_tol);
    EXPECT_NEAR(sc_cst, Ip, m_tol);

    /* -- Change in sign of slope: x, y, z -- */
    // Set up field
    int n = 0;
    init_scalar_slopechange(sc, n, i);
    // Compute interpolated quantities at each face
    get_output_upwind(sc, umac, dt, i, j, k, Im, Ip);
    // Check values
    EXPECT_NEAR(0, Im, m_tol);
    EXPECT_NEAR(0, Ip, m_tol);

    // Set up field
    n = 1;
    init_scalar_slopechange(sc, n, i);
    // Compute interpolated quantities at each face
    get_output_upwind(sc, vmac, dt, i, j, k, Im, Ip);
    // Check values
    EXPECT_NEAR(0, Im, m_tol);
    EXPECT_NEAR(0, Ip, m_tol);

    // Set up field
    n = 2;
    init_scalar_slopechange(sc, n, i);
    // Compute interpolated quantities at each face
    get_output_upwind(sc, wmac, dt, i, j, k, Im, Ip);
    // Check values
    EXPECT_NEAR(0, Im, m_tol);
    EXPECT_NEAR(0, Ip, m_tol);
}

TEST_F(MFluxSchemeTest, minmod)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Parameters to be referenced during test
    amrex::Real adv_vel = 1.2_rt;
    amrex::Real dt = 0.25_rt;

    // Initialize variables to be reused
    amrex::Real Im = 0.0_rt;
    amrex::Real Ip = 0.0_rt;

    // Initialize field variable
    auto& sc = repo.declare_field("scalar", 1, 2);

    // Initialize mac velocity
    repo.declare_face_normal_field({"umac", "vmac", "wmac"}, 1, 1, 1);
    auto& umac = repo.get_field("umac");
    auto& vmac = repo.get_field("vmac");
    auto& wmac = repo.get_field("wmac");
    umac.setVal(adv_vel);
    vmac.setVal(adv_vel);
    wmac.setVal(-adv_vel);

    /* -- Constant field portion of test -- */
    // Set up field
    amrex::Real sc_cst = 1.5_rt;
    init_scalar_uniform(sc, sc_cst);
    // Compute interpolated quantities at each face
    int i = m_nx / 2;
    int j = m_nx / 2;
    int k = m_nx / 2;
    get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
    // Check values
    EXPECT_NEAR(sc_cst, Im, m_tol);
    EXPECT_NEAR(sc_cst, Ip, m_tol);

    /* -- Increasing in slope: x, y, z -- */
    // Values for checking
    auto ir = (amrex::Real)i;
    amrex::Real dx = sc.repo().mesh().Geom(0).CellSizeArray()[0];
    amrex::Real slp =
        (std::pow(ir, 2.0_rt) - std::pow(ir - 1.0_rt, 2.0_rt)) / dx;
    amrex::Real val_p =
        std::pow(ir, 2.0_rt) - (slp * 0.5_rt * (dt * adv_vel - dx));
    amrex::Real val_n =
        std::pow(ir, 2.0_rt) + (slp * 0.5_rt * (dt * adv_vel - dx));
    // Set up field (x)
    init_scalar_increasing(sc, 0);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
    EXPECT_NEAR(i * i, Im, m_tol);
    EXPECT_NEAR(val_p, Ip, m_tol);
    // Set up field (y)
    init_scalar_increasing(sc, 1);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, vmac, dt, i, j, k, 1, Im, Ip);
    EXPECT_NEAR(i * i, Im, m_tol);
    EXPECT_NEAR(val_p, Ip, m_tol);
    // Set up field (z)
    init_scalar_increasing(sc, 2);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, wmac, dt, i, j, k, 2, Im, Ip);
    EXPECT_NEAR(val_n, Im, m_tol);
    EXPECT_NEAR(i * i, Ip, m_tol);

    /* -- Change in sign of slope -- */
    // Set up field (x)
    init_scalar_slopechange(sc, 0, i);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
    // Check values
    EXPECT_NEAR(0, Im, m_tol);
    EXPECT_NEAR(0, Ip, m_tol);
}

TEST_F(MFluxSchemeTest, minmodbdy)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Parameters to be referenced during test
    amrex::Real adv_vel = 1.2_rt;
    amrex::Real dt = 0.25_rt;

    // Initialize variables to be reused
    amrex::Real Im = 0.0_rt;
    amrex::Real Ip = 0.0_rt;

    // Initialize field variable
    auto& sc = repo.declare_field("scalar", 1, 2);

    // Initialize mac velocity
    repo.declare_face_normal_field({"umac", "vmac", "wmac"}, 1, 1, 1);
    auto& umac = repo.get_field("umac");
    auto& vmac = repo.get_field("vmac");
    umac.setVal(-adv_vel);
    vmac.setVal(adv_vel);

    // Initialize BCs
    sc.set_default_fillpatch_bc(sim().time());
    // Set BCs of interest to hoextrap
    sc.bcrec()[0].setLo(0, amrex::BCType::hoextrap);
    sc.bcrec()[0].setHi(1, amrex::BCType::hoextrap);
    // Copy to device to be used
    sc.copy_bc_to_device();

    {
        // Look at lo boundary behavior
        int i = 0;
        int j = m_nx / 2;
        int k = m_nx / 2;
        // Values for checking
        auto ir = (amrex::Real)i;
        amrex::Real dx = sc.repo().mesh().Geom(0).CellSizeArray()[0];
        amrex::Real slp =
            (std::pow(ir, 2.0_rt) - std::pow(ir - 1.0_rt, 2.0_rt)) / dx;
        amrex::Real val_n =
            std::pow(ir, 2.0_rt) + (slp * 0.5_rt * (dt * adv_vel - dx));
        // Set up field
        init_scalar_increasing(sc, 0);
        // Compute interpolated quantities at each face
        get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
        EXPECT_NEAR(val_n, Im, m_tol);
        EXPECT_NEAR(ir * ir, Ip, m_tol);
    }
    {
        // Look at hi boundary behavior
        int i = m_nx / 2;
        int j = m_nx - 1;
        int k = m_nx / 2;
        // Values for checking
        auto ir = static_cast<amrex::Real>(j);
        amrex::Real dx = sc.repo().mesh().Geom(0).CellSizeArray()[0];
        amrex::Real slp =
            (std::pow(ir + 1.0_rt, 2.0_rt) - std::pow(ir, 2.0_rt)) / dx;
        amrex::Real val_p =
            std::pow(ir, 2.0_rt) - (slp * 0.5_rt * (dt * adv_vel - dx));
        // Set up field
        init_scalar_increasing(sc, 1);
        // Compute interpolated quantities at each face
        get_output_minmod(sc, vmac, dt, i, j, k, 1, Im, Ip);
        EXPECT_NEAR(ir * ir, Im, m_tol);
        EXPECT_NEAR(val_p, Ip, m_tol);
    }
}

} // namespace amr_wind_tests
