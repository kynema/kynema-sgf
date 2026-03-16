#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/diagnostics.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_velocity(amr_wind::Field& velocity)
{
    const auto& mesh = velocity.repo().mesh();
    const int nlevels = velocity.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = velocity(lev).arrays();

        amrex::ParallelFor(
            velocity(lev), velocity.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real xc = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real yc = problo[1] + ((j + 0.5_rt) * dx[1]);
                const amrex::Real zc = problo[2] + ((k + 0.5_rt) * dx[2]);

                farrs[nbx](i, j, k, 0) = 1.0_rt - std::pow(xc, 2.0_rt);
                farrs[nbx](i, j, k, 1) = -1.0_rt + std::pow(zc, 2.0_rt);
                farrs[nbx](i, j, k, 2) = 5.0_rt * std::cos(yc);

                if (lev == 0 && nlevels > 1) {
                    // Set base level to large values to detect masking errors
                    farrs[nbx](i, j, k, 0) = 1.0e5_rt;
                    farrs[nbx](i, j, k, 1) = -1.0e5_rt;
                    farrs[nbx](i, j, k, 2) = 1.0e5_rt;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_mac_velocity(
    amr_wind::Field& cc,
    amr_wind::Field& umac,
    amr_wind::Field& vmac,
    amr_wind::Field& wmac)
{
    const auto& mesh = cc.repo().mesh();
    const int nlevels = cc.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& uarrs = umac(lev).arrays();
        const auto& varrs = vmac(lev).arrays();
        const auto& warrs = wmac(lev).arrays();

        amrex::ParallelFor(
            cc(lev), cc.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + (i * dx[0]);
                const amrex::Real yc = problo[1] + ((j + 0.5_rt) * dx[1]);
                const amrex::Real zc = problo[2] + ((k + 0.5_rt) * dx[2]);

                uarrs[nbx](i, j, k) = 1.0_rt - std::pow(x, 2.0_rt);
                varrs[nbx](i, j, k) = -1.0_rt + std::pow(zc, 2.0_rt);
                warrs[nbx](i, j, k) = -3.0_rt * std::cos(yc);

                if (lev == 0 && nlevels > 1) {
                    // Set base level to large values to detect masking errors
                    uarrs[nbx](i, j, k) = 1.0e5_rt;
                    varrs[nbx](i, j, k) = -1.0e5_rt;
                    warrs[nbx](i, j, k) = 1.0e5_rt;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_vof(amr_wind::Field& vof, bool bounded)
{
    const auto& mesh = vof.repo().mesh();
    const int nlevels = vof.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& probhi = mesh.Geom(lev).ProbHiArray();
        const auto& farrs = vof(lev).arrays();

        amrex::ParallelFor(
            vof(lev), vof.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real xc_rel = (i + 0.5_rt) * dx[0];
                const amrex::Real yc_rel = (j + 0.5_rt) * dx[1];
                const amrex::Real zc_rel = (k + 0.5_rt) * dx[2];

                const amrex::Real Lx = probhi[0] - problo[0];
                const amrex::Real Ly = probhi[1] - problo[1];
                const amrex::Real Lz = probhi[2] - problo[2];

                farrs[nbx](i, j, k) =
                    (xc_rel / Lx) + (yc_rel / Ly) + (zc_rel / Lz);

                if (bounded) {
                    farrs[nbx](i, j, k) = amrex::min<amrex::Real>(
                        1.0_rt,
                        amrex::max<amrex::Real>(0.0_rt, farrs[nbx](i, j, k)));
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void modify_vof(amr_wind::Field& vof, amrex::Vector<int> ncell)
{
    const int nlevels = vof.repo().num_active_levels();
    const int nx = ncell[0];
    const int ny = ncell[1];
    const int nz = ncell[2];

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = vof(lev).arrays();

        amrex::ParallelFor(
            vof(lev), vof.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (i == 0 || j == 0 || k == 0) {
                    farrs[nbx](i, j, k) = 0;
                } else if (i == nx - 1 || j == ny - 1 || k == nz - 1) {
                    farrs[nbx](i, j, k) = 1;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

class DiagnosticsTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            pp.addarr("n_cell", m_ncell);
        }
        {
            amrex::ParmParse pp("geometry");

            pp.addarr("prob_lo", m_problo);
            pp.addarr("prob_hi", m_probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    const amrex::Vector<int> m_ncell{{24, 24, 8}};
    const amrex::Vector<amrex::Real> m_problo{{-5.0_rt, -5.0_rt, -2.0_rt}};
    const amrex::Vector<amrex::Real> m_probhi{{5.0_rt, 5.0_rt, 2.0_rt}};
};

TEST_F(DiagnosticsTest, Max_Vel)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& velocity = repo.declare_field("velocity", 3, 0);
    init_velocity(velocity);

    auto cc_results =
        amr_wind::diagnostics::PrintMaxVelLocations(repo, "cell-centered");

    // Check max's and min's, according to profiles
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;
    // max(u)
    EXPECT_NEAR(
        cc_results[0], 1.0_rt - std::pow(0.5_rt * 10.0_rt / 24.0_rt, 2.0_rt),
        tol);
    // min(u)
    EXPECT_NEAR(
        cc_results[4], 1.0_rt - std::pow(11.5_rt * 10.0_rt / 24.0_rt, 2.0_rt),
        tol);
    // max(v)
    EXPECT_NEAR(
        cc_results[8], -1.0_rt + std::pow(3.5_rt * 4.0_rt / 8.0_rt, 2.0_rt),
        tol);
    // min(v)
    EXPECT_NEAR(
        cc_results[12], -1.0_rt + std::pow(0.5_rt * 4.0_rt / 8.0_rt, 2.0_rt),
        tol);
    // max(w)
    EXPECT_NEAR(
        cc_results[16], 5.0_rt * std::cos(0.5_rt * 10.0_rt / 24.0_rt), tol);

    // Check locations (abs due to symmetry)
    EXPECT_NEAR(std::abs(cc_results[1]), 0.5_rt * 10.0_rt / 24.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[5]), 11.5_rt * 10.0_rt / 24.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[11]), 3.5_rt * 4.0_rt / 8.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[15]), 0.5_rt * 4.0_rt / 8.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[18]), 0.5_rt * 10.0_rt / 24.0_rt, tol);
}

TEST_F(DiagnosticsTest, Max_MACvel)
{
    initialize_mesh();
    auto& repo = sim().repo();
    repo.declare_face_normal_field({"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
    auto& umac = repo.get_field("u_mac");
    auto& vmac = repo.get_field("v_mac");
    auto& wmac = repo.get_field("w_mac");
    auto& cc = repo.declare_field("cc", 1, 0);
    init_mac_velocity(cc, umac, vmac, wmac);

    auto fc_results =
        amr_wind::diagnostics::PrintMaxMACVelLocations(repo, "face-centered");

    // Check max's and min's, according to profiles
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;
    // max(umac)
    EXPECT_NEAR(
        fc_results[0], 1.0_rt - std::pow(0.0_rt * 10.0_rt / 24.0_rt, 2.0_rt),
        tol);
    // min(umac)
    EXPECT_NEAR(
        fc_results[4], 1.0_rt - std::pow(12.0_rt * 10.0_rt / 24.0_rt, 2.0_rt),
        tol);
    // max(vmac)
    EXPECT_NEAR(
        fc_results[8], -1.0_rt + std::pow(3.5_rt * 4.0_rt / 8.0_rt, 2.0_rt),
        tol);
    // min(vmac)
    EXPECT_NEAR(
        fc_results[12], -1.0_rt + std::pow(0.5_rt * 4.0_rt / 8.0_rt, 2.0_rt),
        tol);
    // min(wmac)
    EXPECT_NEAR(
        fc_results[20], -3.0_rt * std::cos(0.5_rt * 10.0_rt / 24.0_rt), tol);

    // Check locations
    EXPECT_NEAR(fc_results[1], 0.0_rt * 10.0_rt / 24.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[5]), 12.0_rt * 10.0_rt / 24.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[11]), 3.5_rt * 4.0_rt / 8.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[15]), 0.5_rt * 4.0_rt / 8.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[22]), 0.5_rt * 10.0_rt / 24.0_rt, tol);
}

TEST_F(DiagnosticsTest, Max_Vel_MultiLevel)
{
    populate_parameters();
    {
        amrex::ParmParse pp("amr");
        pp.add("max_level", 1);
    }
    // Create the refinement input file
    // Cover the whole domain for easier testing
    std::stringstream ss;
    ss << "1 // Number of levels" << '\n';
    ss << "1 // Number of boxes at this level" << '\n';
    ss << "-5 -5 -2 5 5 2" << '\n';
    create_mesh_instance<RefineMesh>();
    auto& ref_vec = mesh<RefineMesh>()->refine_criteria_vec();
    ref_vec.emplace_back(std::make_unique<amr_wind::CartBoxRefinement>(sim()));
    auto* box_refine =
        dynamic_cast<amr_wind::CartBoxRefinement*>(ref_vec[0].get());
    box_refine->read_inputs(mesh(), ss);
    initialize_mesh();

    auto& repo = sim().repo();
    auto& velocity = repo.declare_field("velocity", 3, 0);
    init_velocity(velocity);

    auto cc_results =
        amr_wind::diagnostics::PrintMaxVelLocations(repo, "cell-centered");

    // Check max's and min's, according to profiles
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e7_rt;
    // max(u)
    EXPECT_NEAR(
        cc_results[0], 1.0_rt - std::pow(0.5_rt * 10.0_rt / 48.0_rt, 2.0_rt),
        tol);
    // min(u)
    EXPECT_NEAR(
        cc_results[4], 1.0_rt - std::pow(23.5_rt * 10.0_rt / 48.0_rt, 2.0_rt),
        tol);
    // max(v)
    EXPECT_NEAR(
        cc_results[8], -1.0_rt + std::pow(7.5_rt * 4.0_rt / 16.0_rt, 2.0_rt),
        tol);
    // min(v)
    EXPECT_NEAR(
        cc_results[12], -1.0_rt + std::pow(0.5_rt * 4.0_rt / 16.0_rt, 2.0_rt),
        tol);
    // max(w)
    EXPECT_NEAR(
        cc_results[16], 5.0_rt * std::cos(0.5_rt * 10.0_rt / 48.0_rt), tol);

    // Check locations (abs due to symmetry)
    EXPECT_NEAR(std::abs(cc_results[1]), 0.5_rt * 10.0_rt / 48.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[5]), 23.5_rt * 10.0_rt / 48.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[11]), 7.5_rt * 4.0_rt / 16.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[15]), 0.5_rt * 4.0_rt / 16.0_rt, tol);
    EXPECT_NEAR(std::abs(cc_results[18]), 0.5_rt * 10.0_rt / 48.0_rt, tol);
}

TEST_F(DiagnosticsTest, Max_MACvel_MultiLevel)
{
    populate_parameters();
    {
        amrex::ParmParse pp("amr");
        pp.add("max_level", 1);
    }
    // Create the refinement input file
    // Cover the whole domain for easier testing
    std::stringstream ss;
    ss << "1 // Number of levels" << '\n';
    ss << "1 // Number of boxes at this level" << '\n';
    ss << "-5 -5 -2 5 5 2" << '\n';
    create_mesh_instance<RefineMesh>();
    auto& ref_vec = mesh<RefineMesh>()->refine_criteria_vec();
    ref_vec.emplace_back(std::make_unique<amr_wind::CartBoxRefinement>(sim()));
    auto* box_refine =
        dynamic_cast<amr_wind::CartBoxRefinement*>(ref_vec[0].get());
    box_refine->read_inputs(mesh(), ss);
    initialize_mesh();

    auto& repo = sim().repo();
    repo.declare_face_normal_field({"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
    auto& umac = repo.get_field("u_mac");
    auto& vmac = repo.get_field("v_mac");
    auto& wmac = repo.get_field("w_mac");
    auto& cc = repo.declare_field("cc", 1, 0);
    init_mac_velocity(cc, umac, vmac, wmac);

    auto fc_results =
        amr_wind::diagnostics::PrintMaxMACVelLocations(repo, "face-centered");

    // Check max's and min's, according to profiles
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e7_rt;
    // max(umac)
    EXPECT_NEAR(
        fc_results[0], 1.0_rt - std::pow(0.0_rt * 10.0_rt / 48.0_rt, 2.0_rt),
        tol);
    // min(umac)
    EXPECT_NEAR(
        fc_results[4], 1.0_rt - std::pow(24 * 10.0_rt / 48.0_rt, 2.0_rt), tol);
    // max(vmac)
    EXPECT_NEAR(
        fc_results[8], -1.0_rt + std::pow(7.5_rt * 4.0_rt / 16.0_rt, 2.0_rt),
        tol);
    // min(vmac)
    EXPECT_NEAR(
        fc_results[12], -1.0_rt + std::pow(0.5_rt * 4.0_rt / 16.0_rt, 2.0_rt),
        tol);
    // min(wmac)
    EXPECT_NEAR(
        fc_results[20], -3.0_rt * std::cos(0.5_rt * 10.0_rt / 48.0_rt), tol);

    // Check locations
    EXPECT_NEAR(fc_results[1], 0.0_rt * 10.0_rt / 24.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[5]), 24.0_rt * 10.0_rt / 48.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[11]), 7.5_rt * 4.0_rt / 16.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[15]), 0.5_rt * 4.0_rt / 16.0_rt, tol);
    EXPECT_NEAR(std::abs(fc_results[22]), 0.5_rt * 10.0_rt / 48.0_rt, tol);
}

TEST_F(DiagnosticsTest, Field_Extrema)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 1);
    init_vof(vof, false);

    // Get max and min regardless of phase
    amrex::Real fmin{0.}, fmax{0.};
    amr_wind::diagnostics::get_field_extrema(fmax, fmin, vof, 0, 1, 1);

    // Lowest and highest possible values according to init_vof
    const amrex::Real gold_fmin =
        -0.5_rt *
        (1.0_rt / m_ncell[0] + 1.0_rt / m_ncell[1] + 1.0_rt / m_ncell[2]);
    const amrex::Real gold_fmax =
        3.0_rt + (0.5_rt * (1.0_rt / m_ncell[0] + 1.0_rt / m_ncell[1] +
                            1.0_rt / m_ncell[2]));

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;
    EXPECT_NEAR(fmin, gold_fmin, tol);
    EXPECT_NEAR(fmax, gold_fmax, tol);

    // Phase-specific reference values
    const amrex::Real gold_fmin_g = 0.0_rt;
    const amrex::Real gold_fmax_g = 0.0_rt;
    const amrex::Real gold_fmin_l = 1.0_rt;
    const amrex::Real gold_fmax_l = 1.0_rt;

    // Get max and min using masking for phase
    amrex::Real fmin_g{0.0_rt}, fmax_g{0.0_rt}, fmin_l{0.0_rt}, fmax_l{0.0_rt};
    bool found_g = amr_wind::diagnostics::get_field_extrema(
        fmax_g, fmin_g, vof, vof, 0.0_rt, 0, 1, 1);
    bool found_l = amr_wind::diagnostics::get_field_extrema(
        fmax_l, fmin_l, vof, vof, 1.0_rt, 0, 1, 1);

    // Confirm that gas and liquid not found
    EXPECT_FALSE(found_g);
    EXPECT_FALSE(found_l);

    // No gas or liquid in domain, so these values will be bad
    EXPECT_GT(fmin_g, gold_fmin_g);
    EXPECT_LT(fmax_g, gold_fmax_g);
    EXPECT_GT(fmin_l, gold_fmin_l);
    EXPECT_LT(fmax_l, gold_fmax_l);

    // Modify vof to ensure some single-phase cells while remaining unbounded
    modify_vof(vof, m_ncell);
    // Get max and min using masking for phase
    found_g = amr_wind::diagnostics::get_field_extrema(
        fmax_g, fmin_g, vof, vof, 0.0_rt, 0, 1, 1);
    found_l = amr_wind::diagnostics::get_field_extrema(
        fmax_l, fmin_l, vof, vof, 1.0_rt, 0, 1, 1);

    // Phases should be found this time
    EXPECT_TRUE(found_g);
    EXPECT_TRUE(found_l);

    // Extrema should be valid this time
    EXPECT_NEAR(fmin_g, gold_fmin_g, tol);
    EXPECT_NEAR(fmax_g, gold_fmax_g, tol);
    EXPECT_NEAR(fmin_l, gold_fmin_l, tol);
    EXPECT_NEAR(fmax_l, gold_fmax_l, tol);

    // Init vof again with bounds
    init_vof(vof, true);

    // Get max and min with masking on bounded vof
    amrex::Real fmin_g_bounded{0.0_rt}, fmax_g_bounded{0.0_rt},
        fmin_l_bounded{0.0_rt}, fmax_l_bounded{0.0_rt};
    amr_wind::diagnostics::get_field_extrema(
        fmax_g_bounded, fmin_g_bounded, vof, vof, 0.0_rt, 0, 1, 1);
    amr_wind::diagnostics::get_field_extrema(
        fmax_l_bounded, fmin_l_bounded, vof, vof, 1.0_rt, 0, 1, 1);

    EXPECT_NEAR(fmin_g_bounded, gold_fmin_g, tol);
    EXPECT_NEAR(fmax_g_bounded, gold_fmax_g, tol);
    EXPECT_NEAR(fmin_l_bounded, gold_fmin_l, tol);
    EXPECT_NEAR(fmax_l_bounded, gold_fmax_l, tol);
}

} // namespace amr_wind_tests
