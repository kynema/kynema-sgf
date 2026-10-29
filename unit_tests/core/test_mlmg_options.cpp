#include "ks_test_utils/AmrexTest.H"
#include "src/core/MLMGOptions.H"

#include "AMReX_MLPoisson.H"
#include "AMReX_MLABecLaplacian.H"
#include "AMReX_MLNodeLaplacian.H"
#include "AMReX_ParReduce.H"

namespace kynema_sgf_tests {

namespace {

// Domain [0,1] x [0,1] x [0,H] with 16^3 cells: dx/dz = 8
constexpr int ncell = 16;
constexpr amrex::Real height = 0.125_rt;

// max|phi - phi_exact| of the exact discrete solutions on this mesh. The
// operators are diagonal in the Fourier/sine modes of phi, so these follow
// from the discrete eigenvalues: the 7-point cell stencil with second-order
// Dirichlet ghost cells, and the 27-point trilinear nodal stencil with a
// lumped right-hand side (the 1D factors are (2 - 2 cos t) / h^2 and
// (2 + cos t) / 3).
constexpr amrex::Real cell_disc_error = 4.107836e-3_rt;
constexpr amrex::Real nodal_disc_error = 5.500134e-2_rt;

// Solver tolerances and the matching bound on |error - disc_error|. The
// smallest eigenvalue magnitude is about 630 and max|rhs| is about 710, so the
// solver error is about mg_rtol.
#ifdef AMREX_USE_FLOAT
constexpr amrex::Real mg_rtol = 1.0e-5_rt;
constexpr amrex::Real mg_atol = 1.0e-6_rt;
constexpr amrex::Real err_tol = 1.0e-3_rt;
#else
constexpr amrex::Real mg_rtol = 1.0e-10_rt;
constexpr amrex::Real mg_atol = 1.0e-12_rt;
constexpr amrex::Real err_tol = 1.0e-6_rt;
#endif
constexpr int max_iter = 50;

amrex::Geometry make_geom()
{
    const amrex::Box domain(
        amrex::IntVect(0), amrex::IntVect(ncell - 1), amrex::IntVect(0));
    const amrex::RealBox rb({0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, height});
    const amrex::Array<int, AMREX_SPACEDIM> is_periodic{1, 1, 0};
    return amrex::Geometry(domain, rb, amrex::CoordSys::cartesian, is_periodic);
}

amrex::BoxArray make_ba(const amrex::Geometry& geom)
{
    amrex::BoxArray ba(geom.Domain());
    ba.maxSize(ncell / 2);
    return ba;
}

// phi = sin(2 pi x) sin(2 pi y) sin(pi z / H): periodic in x and y, zero at
// z = 0 and z = H. Its Laplacian is -(8 pi^2 + pi^2 / H^2) phi.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
exact_phi(const amrex::Real x, const amrex::Real y, const amrex::Real z)
{
    constexpr amrex::Real pi = amrex::Math::pi<amrex::Real>();
    return std::sin(2.0_rt * pi * x) * std::sin(2.0_rt * pi * y) *
           std::sin(pi * z / height);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
exact_lap(const amrex::Real x, const amrex::Real y, const amrex::Real z)
{
    constexpr amrex::Real pi = amrex::Math::pi<amrex::Real>();
    return -(8.0_rt * pi * pi + pi * pi / (height * height)) *
           exact_phi(x, y, z);
}

// Fill the right-hand side with the exact Laplacian at cell centers
// (offset = 0.5) or nodes (offset = 0)
void fill_rhs(
    amrex::MultiFab& rhs,
    const amrex::Geometry& geom,
    const amrex::Real offset,
    const amrex::Real scale = 1.0_rt)
{
    const auto dx = geom.CellSizeArray();
    const auto& arrs = rhs.arrays();
    amrex::ParallelFor(
        rhs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            arrs[nbx](i, j, k) =
                scale * exact_lap(
                            (i + offset) * dx[0], (j + offset) * dx[1],
                            (k + offset) * dx[2]);
        });
    amrex::Gpu::streamSynchronize();
}

// Max |phi - phi_exact| over the valid region
amrex::Real max_error(
    const amrex::MultiFab& phi,
    const amrex::Geometry& geom,
    const amrex::Real offset,
    const amrex::Real scale = 1.0_rt)
{
    const auto dx = geom.CellSizeArray();
    const auto& arrs = phi.const_arrays();
    amrex::Real err = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMax>{}, amrex::TypeList<amrex::Real>{},
        phi, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Real> {
            const amrex::Real ex =
                scale * exact_phi(
                            (i + offset) * dx[0], (j + offset) * dx[1],
                            (k + offset) * dx[2]);
            return {std::abs(arrs[nbx](i, j, k) - ex)};
        });
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

struct SolveResult
{
    amrex::Real error;
    int num_iters;
    bool used_algmg;
};

// Cell-centered Poisson solve, as used by the MAC projection
SolveResult
cell_solve(const std::string& prefix, const amrex::Real scale = 1.0_rt)
{
    const auto geom = make_geom();
    const auto ba = make_ba(geom);
    const amrex::DistributionMapping dm(ba);

    kynema_sgf::MLMGOptions options(prefix);
    amrex::MLPoisson linop({geom}, {ba}, {dm}, options.lpinfo());
    linop.setMaxOrder(options.max_order);
    linop.setDomainBC(
        {amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
         amrex::LinOpBCType::Dirichlet},
        {amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
         amrex::LinOpBCType::Dirichlet});

    amrex::MultiFab phi(ba, dm, 1, 1);
    amrex::MultiFab rhs(ba, dm, 1, 0);
    phi.setVal(0.0_rt);
    fill_rhs(rhs, geom, 0.5_rt, scale);
    linop.setLevelBC(0, &phi);

    amrex::MLMG mlmg(linop);
    options(mlmg);
    mlmg.solve({&phi}, {&rhs}, options.rel_tol, options.abs_tol);

    return {
        max_error(phi, geom, 0.5_rt, scale), mlmg.getNumIters(),
        mlmg.usedAlgMG()};
}

// Nodal Poisson solve, as used by the nodal projection
SolveResult nodal_solve(const std::string& prefix)
{
    const auto geom = make_geom();
    const auto ba = make_ba(geom);
    const amrex::DistributionMapping dm(ba);
    const auto nba = amrex::convert(ba, amrex::IntVect(1));

    kynema_sgf::MLMGOptions options(prefix);
    amrex::MLNodeLaplacian linop({geom}, {ba}, {dm}, options.lpinfo());
    linop.setDomainBC(
        {amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
         amrex::LinOpBCType::Dirichlet},
        {amrex::LinOpBCType::Periodic, amrex::LinOpBCType::Periodic,
         amrex::LinOpBCType::Dirichlet});
    amrex::MultiFab sigma(ba, dm, 1, 1);
    sigma.setVal(1.0_rt);
    linop.setSigma(0, sigma);

    amrex::MultiFab phi(nba, dm, 1, 1);
    amrex::MultiFab rhs(nba, dm, 1, 0);
    // Dirichlet nodes at z = 0 and z = H keep phi = 0
    phi.setVal(0.0_rt);
    fill_rhs(rhs, geom, 0.0_rt);

    amrex::MLMG mlmg(linop);
    options(mlmg);
    mlmg.solve({&phi}, {&rhs}, options.rel_tol, options.abs_tol);

    return {max_error(phi, geom, 0.0_rt), mlmg.getNumIters(), mlmg.usedAlgMG()};
}

void add_solver_inputs(
    const std::string& prefix, const std::string& multigrid_type)
{
    amrex::ParmParse pp(prefix);
    pp.add("multigrid_type", multigrid_type);
    pp.add("maxiter", max_iter);
    pp.add("mg_rtol", mg_rtol);
    pp.add("mg_atol", mg_atol);
    if (multigrid_type == "hybrid") {
        pp.add("do_fixed_iters", false);
    }
}

amrex::MLMG::MultigridType applied_type(const std::string& prefix)
{
    const auto geom = make_geom();
    const auto ba = make_ba(geom);
    const amrex::DistributionMapping dm(ba);
    kynema_sgf::MLMGOptions options(prefix);
    amrex::MLPoisson linop({geom}, {ba}, {dm}, options.lpinfo());
    amrex::MLMG mlmg(linop);
    options(mlmg);
    return mlmg.getMultigridType();
}

} // namespace

class MLMGOptionsTest : public AmrexTest
{};

TEST_F(MLMGOptionsTest, default_is_geometric)
{
    EXPECT_EQ(applied_type("mac_proj"), amrex::MLMG::MultigridType::geometric);
}

TEST_F(MLMGOptionsTest, multigrid_type_is_applied)
{
    {
        amrex::ParmParse pp("mac_proj");
        pp.add("multigrid_type", std::string("algebraic"));
    }
    {
        amrex::ParmParse pp("nodal_proj");
        pp.add("multigrid_type", std::string("hybrid"));
        pp.add("do_fixed_iters", false);
    }
    {
        amrex::ParmParse pp("diffusion");
        pp.add("multigrid_type", std::string("geometric"));
    }
    EXPECT_EQ(applied_type("mac_proj"), amrex::MLMG::MultigridType::algebraic);
    EXPECT_EQ(applied_type("nodal_proj"), amrex::MLMG::MultigridType::hybrid);
    EXPECT_EQ(applied_type("diffusion"), amrex::MLMG::MultigridType::geometric);
}

TEST_F(MLMGOptionsTest, custom_prefix_overrides_default)
{
    {
        amrex::ParmParse pp("diffusion");
        pp.add("multigrid_type", std::string("algebraic"));
    }
    {
        amrex::ParmParse pp("temperature_diffusion");
        pp.add("multigrid_type", std::string("geometric"));
    }

    const auto geom = make_geom();
    const auto ba = make_ba(geom);
    const amrex::DistributionMapping dm(ba);
    amrex::MLPoisson linop({geom}, {ba}, {dm});
    {
        kynema_sgf::MLMGOptions options("diffusion", "temperature_diffusion");
        amrex::MLMG mlmg(linop);
        options(mlmg);
        EXPECT_EQ(
            mlmg.getMultigridType(), amrex::MLMG::MultigridType::geometric);
    }
    {
        // A custom prefix without the key inherits the default prefix value
        kynema_sgf::MLMGOptions options("diffusion", "tke_diffusion");
        amrex::MLMG mlmg(linop);
        options(mlmg);
        EXPECT_EQ(
            mlmg.getMultigridType(), amrex::MLMG::MultigridType::algebraic);
    }
}

TEST_F(MLMGOptionsTest, bottom_atol_default_under_algmg)
{
    const auto geom = make_geom();
    const auto ba = make_ba(geom);
    const amrex::DistributionMapping dm(ba);
    amrex::MLPoisson linop({geom}, {ba}, {dm});
    {
        amrex::ParmParse pp("mac_proj");
        pp.add("multigrid_type", std::string("algebraic"));
    }
    {
        // Disabled by default with AlgMG
        kynema_sgf::MLMGOptions options("mac_proj");
        amrex::MLMG mlmg(linop);
        options(mlmg);
        EXPECT_LT(mlmg.getBottomToleranceAbs(), 0.0_rt);
    }
    {
        // A value set in the default prefix is kept
        amrex::ParmParse pp("diffusion");
        pp.add("multigrid_type", std::string("algebraic"));
        pp.add("bottom_atol", 1.0e-5_rt);
        kynema_sgf::MLMGOptions options("diffusion", "temperature_diffusion");
        amrex::MLMG mlmg(linop);
        options(mlmg);
        EXPECT_EQ(mlmg.getBottomToleranceAbs(), 1.0e-5_rt);
    }
    {
        // The geometric default is unchanged
        kynema_sgf::MLMGOptions options("nodal_proj");
        amrex::MLMG mlmg(linop);
        options(mlmg);
        EXPECT_GT(mlmg.getBottomToleranceAbs(), 0.0_rt);
    }
}

TEST_F(MLMGOptionsTest, invalid_multigrid_type_aborts)
{
    amrex::ParmParse pp("mac_proj");
    pp.add("multigrid_type", std::string("algebraic_typo"));
    EXPECT_THROW(kynema_sgf::MLMGOptions("mac_proj"), std::runtime_error);
}

TEST_F(MLMGOptionsTest, hybrid_with_fixed_iters_aborts)
{
    // do_fixed_iters defaults to true
    amrex::ParmParse pp("nodal_proj");
    pp.add("multigrid_type", std::string("hybrid"));
    EXPECT_THROW(applied_type("nodal_proj"), std::runtime_error);
}

TEST_F(MLMGOptionsTest, algebraic_with_nsolve_aborts)
{
    amrex::ParmParse pp("nodal_proj");
    pp.add("multigrid_type", std::string("algebraic"));
    pp.add("do_nsolve", true);
    EXPECT_THROW(applied_type("nodal_proj"), std::runtime_error);
}

TEST_F(MLMGOptionsTest, algebraic_with_multicomponent_operator_aborts)
{
    // As in the default velocity diffusion solve
    {
        amrex::ParmParse pp("diffusion");
        pp.add("multigrid_type", std::string("algebraic"));
    }
    const auto geom = make_geom();
    const auto ba = make_ba(geom);
    const amrex::DistributionMapping dm(ba);
    kynema_sgf::MLMGOptions options("diffusion", "velocity_diffusion");
    amrex::MLABecLaplacian linop(
        {geom}, {ba}, {dm}, options.lpinfo(), {}, AMREX_SPACEDIM);
    amrex::MLMG mlmg(linop);
    EXPECT_THROW(options(mlmg), std::runtime_error);
}

TEST_F(MLMGOptionsTest, cell_algebraic_high_aspect_ratio)
{
    add_solver_inputs("mac_proj", "algebraic");
    {
        amrex::ParmParse pp("mac_proj");
        pp.add("max_coarsening_level", 0);
        pp.add("bottom_solver", std::string("amg"));
    }
    const auto res = cell_solve("mac_proj");
    EXPECT_TRUE(res.used_algmg);
    EXPECT_NEAR(res.error, cell_disc_error, err_tol);
    EXPECT_LE(res.num_iters, 3);
}

TEST_F(MLMGOptionsTest, nodal_algebraic_high_aspect_ratio)
{
    add_solver_inputs("nodal_proj", "algebraic");
    {
        amrex::ParmParse pp("nodal_proj");
        pp.add("max_coarsening_level", 0);
        pp.add("bottom_solver", std::string("amg"));
    }
    const auto res = nodal_solve("nodal_proj");
    EXPECT_TRUE(res.used_algmg);
    EXPECT_NEAR(res.error, nodal_disc_error, err_tol);
    EXPECT_LE(res.num_iters, 3);
}

TEST_F(MLMGOptionsTest, cell_geometric_matches_algebraic)
{
    add_solver_inputs("mac_proj", "geometric");
    add_solver_inputs("mac_proj_amg", "algebraic");
    const auto geo = cell_solve("mac_proj");
    const auto alg = cell_solve("mac_proj_amg");
    EXPECT_FALSE(geo.used_algmg);
    EXPECT_TRUE(alg.used_algmg);
    // Geometric multigrid with 2 smoothing sweeps converges slowly at
    // dx/dz = 8 and uses all iterations without meeting mg_rtol (no abort,
    // since do_fixed_iters = true). Its error is within err_tol of the
    // discrete solution after about 36 iterations.
#ifndef AMREX_USE_FLOAT
    EXPECT_EQ(geo.num_iters, max_iter);
#endif
    EXPECT_LT(alg.num_iters, geo.num_iters);
    EXPECT_NEAR(geo.error, cell_disc_error, err_tol);
    EXPECT_NEAR(alg.error, cell_disc_error, err_tol);
    EXPECT_LE(alg.num_iters, 3);
}

#ifndef AMREX_USE_FLOAT
TEST_F(MLMGOptionsTest, cell_algebraic_small_rhs_converges)
{
    // With the default tolerances and a small rhs (max|rhs| about 1e-4, as in
    // the MAC projection after the first step), the default bottom_atol would
    // stop the AlgMG level 0 solve before the MLMG target is met. Without a
    // fixed number of iterations, MLMG then aborts. The single-precision
    // default mg_atol is larger than this rhs, so the test is double only.
    {
        amrex::ParmParse pp("mac_proj");
        pp.add("multigrid_type", std::string("algebraic"));
        pp.add("maxiter", max_iter);
        pp.add("do_fixed_iters", false);
    }
    const amrex::Real scale = 1.0e-7_rt;
    const auto res = cell_solve("mac_proj", scale);
    EXPECT_TRUE(res.used_algmg);
    EXPECT_NEAR(res.error, scale * cell_disc_error, scale * err_tol);
    EXPECT_LE(res.num_iters, 3);
}
#endif

TEST_F(MLMGOptionsTest, cell_hybrid_high_aspect_ratio)
{
    add_solver_inputs("mac_proj", "hybrid");
    const auto res = cell_solve("mac_proj");
#ifndef AMREX_USE_FLOAT
    // Geometric stalls on this mesh before reaching mg_rtol, so hybrid
    // switches to AlgMG
    EXPECT_TRUE(res.used_algmg);
#endif
    EXPECT_NEAR(res.error, cell_disc_error, err_tol);
}

} // namespace kynema_sgf_tests
