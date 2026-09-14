#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/field_boundary_fill/Flather.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {
void initialize_vof(
    kynema_sgf::Field& vof,
    const amrex::Vector<amrex::Geometry>& geom,
    amrex::Real wlev)
{
    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {
        auto& vof_mfab = vof(lev);
        auto vof_arrs = vof_mfab.arrays();
        const auto& zlo = geom[lev].ProbLo(2);
        const auto& dz = geom[lev].CellSize(2);
        amrex::ParallelFor(
            vof_mfab, amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = zlo + (k + 0.5_rt) * dz;

                if (z + 0.5_rt * dz <= wlev) {
                    vof_arrs[nbx](i, j, k) = 1.0_rt;
                } else if (z - 0.5_rt * dz >= wlev) {
                    vof_arrs[nbx](i, j, k) = 0.0_rt;
                } else {
                    vof_arrs[nbx](i, j, k) = (wlev - (z - 0.5_rt * dz)) / dz;
                }
            });
    }
}
} // namespace

class FlatherBoundaryAverageTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{8.0_rt, 8.0_rt, 8.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            amrex::Vector<int> periodic{{0, 0, 0}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("amr");
            const amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 1);
            pp.add("max_grid_size", m_nx);
            pp.add("blocking_factor", 2);
            pp.addarr("n_cell", ncell);
        }

        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "0 0 2 4 8 6" << '\n';

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<kynema_sgf::CartBoxRefinement> box_refine(
            new kynema_sgf::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    const int m_nx{32};
    const amrex::Real m_wlev{4.07_rt};
};

TEST_F(FlatherBoundaryAverageTest, accumulate_boundary_multilevel)
{
    constexpr amrex::Real u0 = 2.0_rt;
    constexpr amrex::Real v0 = 3.0_rt;
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    populate_parameters();
    initialize_mesh();

    auto& repo = mesh().field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 1);
    repo.declare_face_normal_field({"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
    auto& vof = repo.declare_field("vof", 1, 1);

    velocity.setVal(u0, 0, 1, 1);
    velocity.setVal(v0, 1, 1, 1);
    velocity.setVal(0.0_rt, 2, 1, 1);
    initialize_vof(vof, mesh().Geom(), m_wlev);

    kynema_sgf::Flather flather(sim());

    kynema_sgf::MultiLevelVector xlo_uh;
    kynema_sgf::MultiLevelVector xlo_h;
    kynema_sgf::MultiLevelVector xhi_uh;
    kynema_sgf::MultiLevelVector xhi_h;
    kynema_sgf::MultiLevelVector ylo_uh;
    kynema_sgf::MultiLevelVector ylo_h;
    kynema_sgf::MultiLevelVector yhi_uh;
    kynema_sgf::MultiLevelVector yhi_h;

    xlo_uh.resize(1, mesh().Geom());
    xlo_h.resize(1, mesh().Geom());
    xhi_uh.resize(1, mesh().Geom());
    xhi_h.resize(1, mesh().Geom());
    ylo_uh.resize(0, mesh().Geom());
    ylo_h.resize(0, mesh().Geom());
    yhi_uh.resize(0, mesh().Geom());
    yhi_h.resize(0, mesh().Geom());

    const int nlevels = repo.num_active_levels();
    EXPECT_EQ(nlevels, 2);

    // Flow is in positive x and y (outflow at xhi and yhi)
    // - means it is capped at other boundaries
    for (int lev = 0; lev < nlevels; ++lev) {
        flather.accumulate_boundary(
            lev, 0, -1, true, xlo_uh, xlo_h, false,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 0, -1, false, xhi_uh, xhi_h, false,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 1, -1, true, ylo_uh, ylo_h, false,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 1, -1, false, yhi_uh, yhi_h, false,
            kynema_sgf::FieldState::New);

        const auto xhi_idx = xhi_uh.ncells(lev) - 1;
        const auto yhi_idx = yhi_uh.ncells(lev) - 1;

        EXPECT_NEAR(xlo_uh.host_data(lev)[0], 0.0_rt, tol);
        EXPECT_NEAR(xhi_uh.host_data(lev)[xhi_idx], u0 * m_wlev, tol);
        EXPECT_NEAR(ylo_uh.host_data(lev)[0], 0.0_rt, tol);
        EXPECT_NEAR(yhi_uh.host_data(lev)[yhi_idx], v0 * m_wlev, tol);

        EXPECT_NEAR(xlo_h.host_data(lev)[0], m_wlev, tol);
        EXPECT_NEAR(xhi_h.host_data(lev)[xhi_idx], m_wlev, tol);
        EXPECT_NEAR(ylo_h.host_data(lev)[0], m_wlev, tol);
        EXPECT_NEAR(yhi_h.host_data(lev)[yhi_idx], m_wlev, tol);
    }

    velocity.setVal(-u0, 0, 1, 1);
    velocity.setVal(-v0, 1, 1, 1);
    // Flow is now reversed, opposite sides have outflow
    for (int lev = 0; lev < nlevels; ++lev) {
        flather.accumulate_boundary(
            lev, 0, -1, true, xlo_uh, xlo_h, false,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 0, -1, false, xhi_uh, xhi_h, false,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 1, -1, true, ylo_uh, ylo_h, false,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 1, -1, false, yhi_uh, yhi_h, false,
            kynema_sgf::FieldState::New);

        const auto xhi_idx = xhi_uh.ncells(lev) - 1;
        const auto yhi_idx = yhi_uh.ncells(lev) - 1;

        EXPECT_NEAR(xlo_uh.host_data(lev)[0], -u0 * m_wlev, tol);
        EXPECT_NEAR(xhi_uh.host_data(lev)[xhi_idx], 0.0_rt, tol);
        EXPECT_NEAR(ylo_uh.host_data(lev)[0], -v0 * m_wlev, tol);
        EXPECT_NEAR(yhi_uh.host_data(lev)[yhi_idx], 0.0_rt, tol);

        const auto dz = (8.0_rt / (m_nx << lev));

        // Liquid only
        flather.accumulate_boundary(
            lev, 0, 0, true, xlo_uh, xlo_h, false, kynema_sgf::FieldState::New);
        const auto hliq = std::floor(m_wlev / dz) * dz;
        EXPECT_NEAR(xlo_uh.host_data(lev)[0], -u0 * hliq, tol);

        // Mix only
        flather.accumulate_boundary(
            lev, 0, 1, true, xlo_uh, xlo_h, false, kynema_sgf::FieldState::New);
        const auto hmix = m_wlev - hliq;
        EXPECT_NEAR(xlo_uh.host_data(lev)[0], -u0 * hmix, tol);
    }
}

TEST_F(
    FlatherBoundaryAverageTest, accumulate_boundary_multilevel_boundary_cells)
{
    constexpr amrex::Real u0 = 4.0_rt;
    constexpr amrex::Real v0 = 1.5_rt;
    constexpr amrex::Real u0_mac = 11.0_rt;
    constexpr amrex::Real v0_mac = 13.0_rt;
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    populate_parameters();
    initialize_mesh();

    auto& repo = mesh().field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 1);

    const auto mac_fields =
        repo.declare_face_normal_field({"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
    auto& vof = repo.declare_field("vof", 1, 1);

    velocity.setVal(u0, 0, 1, 1);
    velocity.setVal(v0, 1, 1, 1);
    velocity.setVal(0.0_rt, 2, 1, 1);
    mac_fields[0]->setVal(u0_mac, 0, 1, 1);
    mac_fields[1]->setVal(v0_mac, 0, 1, 1);
    mac_fields[2]->setVal(0.0_rt, 0, 1, 1);
    initialize_vof(vof, mesh().Geom(), m_wlev);

    kynema_sgf::Flather flather(sim());

    kynema_sgf::MultiLevelVector xlo_uh;
    kynema_sgf::MultiLevelVector xlo_h;
    kynema_sgf::MultiLevelVector xhi_uh;
    kynema_sgf::MultiLevelVector xhi_h;
    kynema_sgf::MultiLevelVector ylo_uh;
    kynema_sgf::MultiLevelVector ylo_h;
    kynema_sgf::MultiLevelVector yhi_uh;
    kynema_sgf::MultiLevelVector yhi_h;

    xlo_uh.resize(1, mesh().Geom());
    xlo_h.resize(1, mesh().Geom());
    xhi_uh.resize(1, mesh().Geom());
    xhi_h.resize(1, mesh().Geom());
    ylo_uh.resize(0, mesh().Geom());
    ylo_h.resize(0, mesh().Geom());
    yhi_uh.resize(0, mesh().Geom());
    yhi_h.resize(0, mesh().Geom());

    const int nlevels = repo.num_active_levels();
    EXPECT_EQ(nlevels, 2);

    for (int lev = 0; lev < nlevels; ++lev) {
        flather.accumulate_boundary(
            lev, 0, -1, true, xlo_uh, xlo_h, true, kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 0, -1, false, xhi_uh, xhi_h, true,
            kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 1, -1, true, ylo_uh, ylo_h, true, kynema_sgf::FieldState::New);
        flather.accumulate_boundary(
            lev, 1, -1, false, yhi_uh, yhi_h, true,
            kynema_sgf::FieldState::New);

        EXPECT_NEAR(xlo_uh.host_data(lev)[0], u0 * m_wlev, tol);
        EXPECT_NEAR(xhi_uh.host_data(lev)[0], u0 * m_wlev, tol);
        EXPECT_NEAR(ylo_uh.host_data(lev)[0], v0 * m_wlev, tol);
        EXPECT_NEAR(yhi_uh.host_data(lev)[0], v0 * m_wlev, tol);

        EXPECT_NEAR(xlo_h.host_data(lev)[0], m_wlev, tol);
        EXPECT_NEAR(xhi_h.host_data(lev)[0], m_wlev, tol);
        EXPECT_NEAR(ylo_h.host_data(lev)[0], m_wlev, tol);
        EXPECT_NEAR(yhi_h.host_data(lev)[0], m_wlev, tol);

        // MAC velocity
        flather.accumulate_boundary(
            lev, 0, -1, true, xlo_uh, xlo_h, true, kynema_sgf::FieldState::New,
            true);
        flather.accumulate_boundary(
            lev, 1, -1, true, ylo_uh, ylo_h, true, kynema_sgf::FieldState::New,
            true);

        EXPECT_NEAR(xlo_uh.host_data(lev)[0], u0_mac * m_wlev, tol);
        EXPECT_NEAR(ylo_uh.host_data(lev)[0], v0_mac * m_wlev, tol);
    }
}

class FlatherThreeLevelTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{8.0_rt, 8.0_rt, 8.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            amrex::Vector<int> periodic{{0, 0, 0}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("amr");
            const amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 2);
            pp.add("max_grid_size", m_nx);
            pp.add("blocking_factor", 2);
            pp.addarr("n_cell", ncell);
        }

        // Two levels of refinement, nested so that the finest level is two
        // refinement ratios removed from the coarsest
        std::stringstream ss;
        ss << "2 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "0 0 2 4 8 6" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "0 0 3 2 8 5" << '\n';

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<kynema_sgf::CartBoxRefinement> box_refine(
            new kynema_sgf::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    const int m_nx{32};
    const amrex::Real m_wlev{4.07_rt};
};

TEST_F(FlatherThreeLevelTest, accumulate_boundary_index_translation)
{
    constexpr amrex::Real u0 = 2.0_rt;
    constexpr amrex::Real v0 = 3.0_rt;
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    populate_parameters();
    initialize_mesh();

    auto& repo = mesh().field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 1);
    repo.declare_face_normal_field({"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
    auto& vof = repo.declare_field("vof", 1, 1);

    // Negative velocity makes xlo and ylo the outflow boundaries, so the
    // interior sums are not clipped to zero
    velocity.setVal(-u0, 0, 1, 1);
    velocity.setVal(-v0, 1, 1, 1);
    velocity.setVal(0.0_rt, 2, 1, 1);
    initialize_vof(vof, mesh().Geom(), m_wlev);

    kynema_sgf::Flather flather(sim());

    const int nlevels = repo.num_active_levels();
    ASSERT_EQ(nlevels, 3);
    const int finest = nlevels - 1;

    kynema_sgf::MultiLevelVector xlo_uh;
    kynema_sgf::MultiLevelVector xlo_h;
    kynema_sgf::MultiLevelVector ylo_uh;
    kynema_sgf::MultiLevelVector ylo_h;

    xlo_uh.resize(1, mesh().Geom());
    xlo_h.resize(1, mesh().Geom());
    ylo_uh.resize(0, mesh().Geom());
    ylo_h.resize(0, mesh().Geom());

    // Accumulating on the finest level draws contributions from levels that
    // are one and two refinement ratios coarser
    flather.accumulate_boundary(
        finest, 0, -1, true, xlo_uh, xlo_h, false, kynema_sgf::FieldState::New);
    flather.accumulate_boundary(
        finest, 1, -1, true, ylo_uh, ylo_h, false, kynema_sgf::FieldState::New);

    ASSERT_EQ(xlo_h.ncells(finest), m_nx * 4);
    ASSERT_EQ(ylo_h.ncells(finest), m_nx * 4);

    // Each column is covered by exactly one level, so a coarse index must map
    // onto a contiguous range of fine indices with no gaps and no overlap.
    // A gap leaves a column short and an overlap double counts it.
    for (int n = 0; n < xlo_h.ncells(finest); ++n) {
        EXPECT_NEAR(xlo_h.host_data(finest)[n], m_wlev, tol)
            << "xlo liquid height, column " << n;
        EXPECT_NEAR(xlo_uh.host_data(finest)[n], -u0 * m_wlev, tol)
            << "xlo depth-integrated velocity, column " << n;
    }
    for (int n = 0; n < ylo_h.ncells(finest); ++n) {
        EXPECT_NEAR(ylo_h.host_data(finest)[n], m_wlev, tol)
            << "ylo liquid height, column " << n;
        EXPECT_NEAR(ylo_uh.host_data(finest)[n], -v0 * m_wlev, tol)
            << "ylo depth-integrated velocity, column " << n;
    }
}

} // namespace kynema_sgf_tests