#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ForestCanopy.H"
#include "src/core/field_ops.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"
#include <cmath>

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

// value below ztop, zero above
void set_below(
    kynema_sgf::Field& fld, const amrex::Real value, const amrex::Real ztop)
{
    const auto& mesh = fld.repo().mesh();
    for (int lev = 0; lev < fld.repo().num_active_levels(); ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                farrs[nbx](i, j, k) = (z < ztop) ? value : 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

class ForestCanopyTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 8}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{8.0_rt, 8.0_rt, 8.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
    }

    // Canopy of forest drag fd below 4 m, velocity (3, 4, 0), uniform tke
    void setup_fields(const amrex::Real fd, const amrex::Real tke_value)
    {
        populate_parameters();
        initialize_mesh();
        auto& repo = sim().repo();
        auto& forest = repo.declare_field("forest_drag", 1, 1, 1);
        auto& vel = repo.declare_field("velocity", 3, 1, 1);
        auto& tke = repo.declare_field("tke", 1, 1, 1);
        auto& src = repo.declare_field("tke_src_term", 1, 0, 1);
        set_below(forest, fd, 4.0_rt);
        vel.setVal(3.0_rt, 0, 1, 1);
        vel.setVal(4.0_rt, 1, 1, 1);
        vel.setVal(0.0_rt, 2, 1, 1);
        tke.setVal(tke_value);
        src.setVal(0.0_rt);
    }

    void add_source(
        const amrex::Real beta_p,
        const amrex::Real beta_d,
        const amrex::Real dt)
    {
        auto& repo = sim().repo();
        kynema_sgf::forestdrag::add_canopy_tke_source(
            repo.get_field("tke_src_term")(0), repo.get_field("forest_drag")(0),
            repo.get_field("velocity")(0), repo.get_field("tke")(0), beta_p,
            beta_d, dt);
    }
};

// Defaults: off, beta_p = 1, beta_d = 4, alpha = 0.04
TEST_F(ForestCanopyTest, defaults)
{
    const auto canopy = kynema_sgf::forestdrag::parse_canopy_turbulence();
    EXPECT_FALSE(canopy.m_enabled);
    EXPECT_EQ(canopy.m_beta_p, 1.0_rt);
    EXPECT_EQ(canopy.m_beta_d, 4.0_rt);
    EXPECT_EQ(canopy.m_length_alpha, 0.04_rt);
}

TEST_F(ForestCanopyTest, inputs)
{
    amrex::ParmParse pp("ForestDrag");
    pp.add("canopy_tke", true);
    pp.add("canopy_beta_p", 0.8_rt);
    pp.add("canopy_beta_d", 5.1_rt);
    pp.add("canopy_length_alpha", 0.1_rt);
    const auto canopy = kynema_sgf::forestdrag::parse_canopy_turbulence();
    EXPECT_TRUE(canopy.m_enabled);
    EXPECT_EQ(canopy.m_beta_p, 0.8_rt);
    EXPECT_EQ(canopy.m_beta_d, 5.1_rt);
    EXPECT_EQ(canopy.m_length_alpha, 0.1_rt);
}

TEST_F(ForestCanopyTest, negative_coefficient_aborts)
{
    amrex::ParmParse pp("ForestDrag");
    pp.add("canopy_beta_d", -1.0_rt);
    EXPECT_THROW(
        kynema_sgf::forestdrag::parse_canopy_turbulence(), amrex::RuntimeError);
}

// S_k = fd (beta_p |U|^3 - beta_d |U| k), with the sink integrated over dt
TEST_F(ForestCanopyTest, tke_source_in_canopy_cell)
{
    const amrex::Real fd = 0.05_rt;
    const amrex::Real tke = 0.5_rt;
    const amrex::Real beta_p = 1.5_rt;
    const amrex::Real beta_d = 3.0_rt;
    const amrex::Real dt = 0.1_rt;
    setup_fields(fd, tke);
    add_source(beta_p, beta_d, dt);

    const amrex::Real ws = 5.0_rt;
    const amrex::Real expected =
        (beta_p * fd * ws * ws * ws) +
        (tke * std::expm1(-beta_d * fd * ws * dt) / dt);
    const auto& src = sim().repo().get_field("tke_src_term");
    // Values of order 10
    const amrex::Real tol = 10.0_rt * kynema_sgf::constants::TIGHT_TOL;
    EXPECT_NEAR(utils::field_probe(src, 0, 3, 3, 0), expected, tol);
    EXPECT_NEAR(utils::field_probe(src, 0, 5, 2, 3), expected, tol);
    // Small beta_d fd |U| dt: close to the explicit value 9.375 - 0.375
    EXPECT_NEAR(expected, 9.0_rt, 0.03_rt);

    // Above the canopy the source is untouched
    EXPECT_EQ(utils::field_probe(src, 0, 3, 3, 4), 0.0_rt);
    EXPECT_EQ(utils::field_probe(src, 0, 3, 3, 7), 0.0_rt);
}

// A dense canopy with a large time step removes at most k: the explicit sink
// -beta_d fd |U| k would make k negative
TEST_F(ForestCanopyTest, stable_sink)
{
    const amrex::Real fd = 10.0_rt;
    const amrex::Real tke = 0.5_rt;
    const amrex::Real beta_d = 4.0_rt;
    const amrex::Real dt = 1.0_rt;
    setup_fields(fd, tke);
    add_source(0.0_rt, beta_d, dt);

    const auto& src = sim().repo().get_field("tke_src_term");
    const amrex::Real sk = utils::field_probe(src, 0, 3, 3, 0);
    const amrex::Real k_new = tke + (dt * sk);
    const amrex::Real explicit_k_new = tke - (dt * beta_d * fd * 5.0_rt * tke);
    EXPECT_LT(explicit_k_new, 0.0_rt);
    EXPECT_GE(k_new, 0.0_rt);
    EXPECT_NEAR(
        k_new, tke * std::exp(-beta_d * fd * 5.0_rt * dt),
        kynema_sgf::constants::TIGHT_TOL);
    EXPECT_GE(sk, -tke / dt);
}

// Small time step: the integrated sink tends to -beta_d fd |U| k
TEST_F(ForestCanopyTest, sink_small_time_step)
{
    const amrex::Real fd = 0.05_rt;
    const amrex::Real tke = 0.5_rt;
    const amrex::Real dt = 1.0e-6_rt;
    setup_fields(fd, tke);
    add_source(0.0_rt, 4.0_rt, dt);
    const auto& src = sim().repo().get_field("tke_src_term");
    const amrex::Real explicit_sink = -4.0_rt * fd * 5.0_rt * tke;
    EXPECT_NEAR(
        utils::field_probe(src, 0, 3, 3, 0), explicit_sink,
        1.0e-5_rt * std::abs(explicit_sink));
}

// Length scale limit l <= alpha / fd, only where fd > 0: l = 10 m, with
// mu_t, shear and buoyancy production rescaled by the same factor
class ForestCanopyLengthTest : public ForestCanopyTest
{
protected:
    void setup_length_fields(const amrex::Real fd)
    {
        setup_fields(fd, 0.5_rt);
        auto& repo = sim().repo();
        repo.declare_field("turb_lscale", 1, 0, 1).setVal(10.0_rt);
        repo.declare_field("mu_turb", 1, 0, 1).setVal(2.0_rt);
        repo.declare_field("shear_prod", 1, 0, 1).setVal(3.0_rt);
        repo.declare_field("buoy_prod", 1, 0, 1).setVal(-1.0_rt);
    }

    void limit(const amrex::Real alpha)
    {
        auto& repo = sim().repo();
        kynema_sgf::forestdrag::limit_canopy_length_scale(
            repo.get_field("forest_drag")(0), alpha,
            repo.get_field("turb_lscale")(0), repo.get_field("mu_turb")(0),
            repo.get_field("shear_prod")(0), repo.get_field("buoy_prod")(0));
    }

    void expect_cell(
        const int k,
        const amrex::Real lscale,
        const amrex::Real ratio,
        const amrex::Real tol)
    {
        auto& repo = sim().repo();
        EXPECT_NEAR(
            utils::field_probe(repo.get_field("turb_lscale"), 0, 3, 3, k),
            lscale, tol);
        EXPECT_NEAR(
            utils::field_probe(repo.get_field("mu_turb"), 0, 3, 3, k),
            2.0_rt * ratio, tol);
        EXPECT_NEAR(
            utils::field_probe(repo.get_field("shear_prod"), 0, 3, 3, k),
            3.0_rt * ratio, tol);
        EXPECT_NEAR(
            utils::field_probe(repo.get_field("buoy_prod"), 0, 3, 3, k),
            -1.0_rt * ratio, tol);
    }
};

TEST_F(ForestCanopyLengthTest, limit_inside_canopy_only)
{
    const amrex::Real fd = 0.02_rt;
    const amrex::Real alpha = 0.04_rt;
    setup_length_fields(fd);
    limit(alpha);

    // alpha / fd = 2 m < 10 m in the canopy (z < 4 m)
    const amrex::Real tol = 10.0_rt * kynema_sgf::constants::TIGHT_TOL;
    const amrex::Real lcap = alpha / fd;
    expect_cell(0, lcap, lcap / 10.0_rt, tol);
    expect_cell(3, lcap, lcap / 10.0_rt, tol);
    // No forest drag: bit-identical
    expect_cell(4, 10.0_rt, 1.0_rt, 0.0_rt);
    expect_cell(7, 10.0_rt, 1.0_rt, 0.0_rt);
}

// A sparse canopy, alpha / fd = 20 m > 10 m, leaves the length scale
TEST_F(ForestCanopyLengthTest, limit_above_length_scale)
{
    setup_length_fields(0.002_rt);
    limit(0.04_rt);
    expect_cell(0, 10.0_rt, 1.0_rt, 0.0_rt);
    expect_cell(3, 10.0_rt, 1.0_rt, 0.0_rt);
}

// A vanishing drag must not overflow alpha / fd
TEST_F(ForestCanopyLengthTest, limit_tiny_drag)
{
    setup_length_fields(1.0e-30_rt);
    limit(0.04_rt);
    expect_cell(0, 10.0_rt, 1.0_rt, 0.0_rt);
}

} // namespace kynema_sgf_tests
