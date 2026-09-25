#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/equation_systems/SchemeTraits.H"
#include "src/equation_systems/vof/vof.H"
#include "src/equation_systems/vof/vof_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

//! 4 x 4 blanked column block within a uniform vof field
void init_blanked_block(
    kynema_sgf::Field& vof,
    kynema_sgf::IntField& blanking,
    const amrex::Real vof_fluid)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& vof_arr = vof(lev).array(mfi);
        const auto& blank_arr = blanking(lev).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const bool in_terrain =
                (i >= 2) && (i <= 5) && (j >= 2) && (j <= 5);
            blank_arr(i, j, k) = in_terrain ? 1 : 0;
            vof_arr(i, j, k) = in_terrain ? 0.0_rt : vof_fluid;
        });
    });
    vof.fillpatch(0.0_rt);
}

//! Single blanked plane at i = 4, with a vof field that varies in x and z
void init_blanked_plane(kynema_sgf::Field& vof, kynema_sgf::IntField& blanking)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& vof_arr = vof(lev).array(mfi);
        const auto& blank_arr = blanking(lev).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const bool in_terrain = (i == 4);
            blank_arr(i, j, k) = in_terrain ? 1 : 0;
            vof_arr(i, j, k) =
                in_terrain ? 0.0_rt : (0.1_rt * i) + (0.01_rt * k);
        });
    });
    vof.fillpatch(0.0_rt);
}

void init_diagonal_case(
    kynema_sgf::Field& vof, kynema_sgf::IntField& blanking)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& vof_arr = vof(lev).array(mfi);
        const auto& blank_arr = blanking(lev).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const bool is_center = (i == 4) && (j == 4);
            const bool is_cardinal =
                ((i == 4) && (j == 3 || j == 5)) ||
                ((j == 4) && (i == 3 || i == 5));
            const bool is_diagonal =
                (i == 3 || i == 5) && (j == 3 || j == 5);
            blank_arr(i, j, k) = (is_center || is_cardinal) ? 1 : 0;
            vof_arr(i, j, k) = is_diagonal ? (0.1_rt * i + 0.01_rt * j)
                                          : 0.0_rt;
        });
    });
    vof.fillpatch(0.0_rt);
}

void init_boundary_source_case(
    kynema_sgf::Field& vof, kynema_sgf::IntField& blanking)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& vof_arr = vof(lev).array(mfi);
        const auto& blank_arr = blanking(lev).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const bool is_target = (i == 1) && (j == 4);
            const bool is_boundary = (i == 0) || (i == 7);
            blank_arr(i, j, k) = is_target ? 1 : 0;
            vof_arr(i, j, k) = is_boundary ? 1.0_rt : 0.2_rt;
        });
    });
    vof.fillpatch(0.0_rt);
}

void get_error_uniform(
    kynema_sgf::ScratchField& err_fld,
    kynema_sgf::Field& vof,
    const amrex::Real vof_fluid)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& err_arr = err_fld(lev).array(mfi);
        const auto& vof_arr = vof(lev).const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            err_arr(i, j, k) = amrex::Math::abs(vof_arr(i, j, k) - vof_fluid);
        });
    });
}

void get_error_blanked_plane(
    kynema_sgf::ScratchField& err_fld, kynema_sgf::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& err_arr = err_fld(lev).array(mfi);
        const auto& vof_arr = vof(lev).const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Average of the i = 3 and i = 5 neighbors within the same plane
            const amrex::Real expected = (i == 4)
                                             ? (0.4_rt + (0.01_rt * k))
                                             : ((0.1_rt * i) + (0.01_rt * k));
            err_arr(i, j, k) = amrex::Math::abs(vof_arr(i, j, k) - expected);
        });
    });
}

void get_error_diagonal_case(
    kynema_sgf::ScratchField& err_fld, kynema_sgf::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& err_arr = err_fld(lev).array(mfi);
        const auto& vof_arr = vof(lev).const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const bool is_center = (i == 4) && (j == 4);
            const bool is_diagonal =
                (i == 3 || i == 5) && (j == 3 || j == 5);
            const amrex::Real expected =
                is_center   ? (0.1_rt * 4.0_rt + 0.01_rt * 4.0_rt)
                : is_diagonal ? (0.1_rt * i + 0.01_rt * j)
                              : 0.0_rt;
            err_arr(i, j, k) =
                (is_center || is_diagonal)
                    ? amrex::Math::abs(vof_arr(i, j, k) - expected)
                    : 0.0_rt;
        });
    });
}

void get_error_boundary_source_case(
    kynema_sgf::ScratchField& err_fld, kynema_sgf::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.validbox();
        const auto& err_arr = err_fld(lev).array(mfi);
        const auto& vof_arr = vof(lev).const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real expected =
                ((i == 1) && (j == 4)) ? 0.5_rt
                                       : (((i == 0) || (i == 7)) ? 1.0_rt
                                                                  : 0.2_rt);
            err_arr(i, j, k) =
                amrex::Math::abs(vof_arr(i, j, k) - expected);
        });
    });
}

class VOFTerrainFillTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            // Multiple boxes to exercise the ghost exchange between sweeps
            amrex::ParmParse pp("amr");
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
        }
        {
            amrex::ParmParse pp("MultiPhase");
            pp.add("density_fluid1", 1000.0_rt);
            pp.add("density_fluid2", 1.0_rt);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"MultiPhase"};
            pp.addarr("physics", physics);
            pp.add("use_godunov", 1);
        }
    }

    void setup_sim()
    {
        initialize_mesh();

        auto& pde_mgr = sim().pde_manager();
        auto& mom_eqn = pde_mgr.register_icns();
        mom_eqn.initialize();

        // Registers the VOF equation and the associated fields
        sim().init_physics();

        sim().repo().declare_int_field("terrain_blank", 1, 1, 1);
    }

    kynema_sgf::pde::PDEFields& vof_fields()
    {
        return sim()
            .pde_manager()(
                kynema_sgf::pde::VOF::pde_name() + "-" +
                kynema_sgf::fvm::Godunov::scheme_name())
            .fields();
    }
};

constexpr amrex::Real tol =
    std::numeric_limits<amrex::Real>::epsilon() * 1.0e2_rt;

} // namespace

//! Blanked cells several layers deep must all recover a uniform vof value
TEST_F(VOFTerrainFillTest, uniform_vof_block)
{
    setup_sim();

    auto& repo = sim().repo();
    auto& vof = repo.get_field("vof");
    auto& blanking = repo.get_int_field("terrain_blank");
    constexpr amrex::Real vof_fluid = 0.7_rt;

    // The innermost blanked cells have no unblanked lateral neighbors and are
    // only reached by successive sweeps
    init_blanked_block(vof, blanking, vof_fluid);

    kynema_sgf::pde::PostSolveOp<kynema_sgf::pde::VOF> post_solve(
        sim(), vof_fields());
    post_solve.extrapolate_vof_into_terrain();

    auto error_ptr = repo.create_scratch_field(1, 0);
    auto& error_fld = *error_ptr;
    get_error_uniform(error_fld, vof, vof_fluid);

    EXPECT_NEAR(error_fld(0).max(0), 0.0_rt, tol);
}

//! Blanked cells average their lateral neighbors only, not their z neighbors
TEST_F(VOFTerrainFillTest, lateral_average_only)
{
    setup_sim();

    auto& repo = sim().repo();
    auto& vof = repo.get_field("vof");
    auto& blanking = repo.get_int_field("terrain_blank");

    init_blanked_plane(vof, blanking);

    kynema_sgf::pde::PostSolveOp<kynema_sgf::pde::VOF> post_solve(
        sim(), vof_fields());
    post_solve.extrapolate_vof_into_terrain();

    auto error_ptr = repo.create_scratch_field(1, 0);
    auto& error_fld = *error_ptr;
    get_error_blanked_plane(error_fld, vof);

    EXPECT_NEAR(error_fld(0).max(0), 0.0_rt, tol);
}

TEST_F(VOFTerrainFillTest, lateral_diagonal_average)
{
    setup_sim();

    auto& repo = sim().repo();
    auto& vof = repo.get_field("vof");
    auto& blanking = repo.get_int_field("terrain_blank");

    init_diagonal_case(vof, blanking);

    kynema_sgf::pde::PostSolveOp<kynema_sgf::pde::VOF> post_solve(
        sim(), vof_fields());
    post_solve.extrapolate_vof_into_terrain();

    auto error_ptr = repo.create_scratch_field(1, 0);
    auto& error_fld = *error_ptr;
    get_error_diagonal_case(error_fld, vof);

    EXPECT_NEAR(error_fld(0).max(0), 0.0_rt, tol);
}

TEST_F(VOFTerrainFillTest, includes_valid_domain_boundary_sources)
{
    setup_sim();

    auto& repo = sim().repo();
    auto& vof = repo.get_field("vof");
    auto& blanking = repo.get_int_field("terrain_blank");

    init_boundary_source_case(vof, blanking);

    kynema_sgf::pde::PostSolveOp<kynema_sgf::pde::VOF> post_solve(
        sim(), vof_fields());
    post_solve.extrapolate_vof_into_terrain();

    auto error_ptr = repo.create_scratch_field(1, 0);
    auto& error_fld = *error_ptr;
    get_error_boundary_source_case(error_fld, vof);

    EXPECT_NEAR(error_fld(0).max(0), 0.0_rt, tol);
}

//! Cells that cannot be reached from any unblanked cell are left alone
TEST_F(VOFTerrainFillTest, fully_blanked_domain)
{
    setup_sim();

    auto& repo = sim().repo();
    auto& vof = repo.get_field("vof");
    auto& blanking = repo.get_int_field("terrain_blank");
    constexpr amrex::Real vof_init = 0.25_rt;

    blanking.setVal(1);
    vof.setVal(vof_init);

    kynema_sgf::pde::PostSolveOp<kynema_sgf::pde::VOF> post_solve(
        sim(), vof_fields());
    post_solve.extrapolate_vof_into_terrain();

    EXPECT_NEAR(vof(0).max(0), vof_init, tol);
    EXPECT_NEAR(vof(0).min(0), vof_init, tol);
}

} // namespace kynema_sgf_tests
