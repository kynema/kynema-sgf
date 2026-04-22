#include <algorithm>

#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/overset/overset_ops_routines.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class VOFOversetOps : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 8}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{1.0_rt, 1.0_rt, 1.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace {

void init_vof_only(kynema_sgf::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 2), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k > 3) {
                    vof_arr(i, j, k) = 0.0_rt;
                } else if (k < 3) {
                    vof_arr(i, j, k) = 1.0_rt;
                } else {
                    vof_arr(i, j, k) = 0.5_rt;
                }
            });
    });
}

void init_iblank_node(kynema_sgf::IntField& iblank)
{
    run_algorithm(iblank, [&](const int lev, const amrex::MFIter& mfi) {
        auto ibl_n = iblank(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 2), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                ibl_n(i, j, k) = 1;
                bool in_overset_x = (i >= 1 && i <= 6);
                bool in_overset_y = (j >= 1 && j <= 6);
                bool in_overset_z = (k >= 1 && k <= 7);
                bool in_block_x = (i >= 3 && i <= 4);
                bool in_block_y = (j >= 3 && j <= 4);
                bool in_block_z = (k >= 3 && k <= 6);
                if (in_overset_x && in_overset_y && in_overset_z) {
                    ibl_n(i, j, k) = -1;
                    if (in_block_x && in_block_y && in_block_z) {
                        ibl_n(i, j, k) = 0;
                    }
                }
            });
    });
}

void init_iblank_cell(kynema_sgf::IntField& iblank)
{
    run_algorithm(iblank, [&](const int lev, const amrex::MFIter& mfi) {
        auto ibl_c = iblank(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(
            grow(bx, 2), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                ibl_c(i, j, k) = 1;
                bool in_overset_x = (i >= 1 && i <= 6);
                bool in_overset_y = (j >= 1 && j <= 6);
                bool in_overset_z = (k >= 1 && k <= 6);
                bool in_block_x = (i >= 3 && i <= 4);
                bool in_block_y = (j >= 3 && j <= 4);
                bool in_block_z = (k >= 3 && k <= 4);
                if (in_overset_x && in_overset_y && in_overset_z) {
                    ibl_c(i, j, k) = -1;
                    if (in_block_x && in_block_y && in_block_z) {
                        ibl_c(i, j, k) = 0;
                    }
                }
            });
    });
}

amrex::Real check_iblank_node_impl(kynema_sgf::IntField& mask_node)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < mask_node.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            mask_node(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<int const> const& i_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    bool in_overset_x = (i >= 1 && i <= 6);
                    bool in_overset_y = (j >= 1 && j <= 6);
                    bool in_overset_z = (k >= 1 && k <= 7);
                    bool near_block_x = (i >= 2 && i <= 5);
                    bool near_block_y = (j >= 2 && j <= 5);
                    bool near_block_z = (k >= 2 && k <= 7);
                    bool in_interface_band = (k >= 2 && k <= 5);
                    if (near_block_x && near_block_y && near_block_z) {
                        error += std::abs(i_arr(i, j, k) - 0);
                    } else {
                        if (!in_interface_band && in_overset_x &&
                            in_overset_y && in_overset_z) {
                            error += std::abs(i_arr(i, j, k) - 0);
                        } else {
                            error += std::abs(i_arr(i, j, k) - 1);
                        }
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real check_iblank_cell_impl(kynema_sgf::IntField& mask_cell)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < mask_cell.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            mask_cell(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<int const> const& i_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    bool in_overset_x = (i >= 1 && i <= 6);
                    bool in_overset_y = (j >= 1 && j <= 6);
                    bool in_overset_z = (k >= 1 && k <= 6);
                    bool near_block_x = (i >= 2 && i <= 5);
                    bool near_block_y = (j >= 2 && j <= 5);
                    bool near_block_z = (k >= 2 && k <= 5);
                    bool in_interface_band = (k >= 2 && k <= 4);
                    if (near_block_x && near_block_y && near_block_z) {
                        error += std::abs(i_arr(i, j, k) - 0);
                    } else {
                        if (!in_interface_band && in_overset_x &&
                            in_overset_y && in_overset_z) {
                            error += std::abs(i_arr(i, j, k) - 0);
                        } else {
                            error += std::abs(i_arr(i, j, k) - 1);
                        }
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real check_iblank_cell_default_impl(kynema_sgf::IntField& mask_cell)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < mask_cell.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            mask_cell(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<int const> const& i_arr) -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    bool in_overset_x = (i >= 1 && i <= 6);
                    bool in_overset_y = (j >= 1 && j <= 6);
                    bool in_overset_z = (k >= 1 && k <= 6);
                    if (in_overset_x && in_overset_y && in_overset_z) {
                        error += std::abs(i_arr(i, j, k) - 0);
                    } else {
                        error += std::abs(i_arr(i, j, k) - 1);
                    }
                });

                return error;
            });
    }
    return error_total;
}
} // namespace

TEST_F(VOFOversetOps, projection_masks)
{
    populate_parameters();
    initialize_mesh();

    auto& repo = sim().repo();
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", 1, nghost);
    auto& iblank_node = repo.declare_int_field(
        "iblank_node", 1, nghost, 1, kynema_sgf::FieldLoc::NODE);
    auto& mask_node = repo.declare_int_field(
        "mask_node", 1, nghost, 1, kynema_sgf::FieldLoc::NODE);
    auto& iblank_cell = repo.declare_int_field("iblank_cell", 1, nghost);
    auto& mask_cell = repo.declare_int_field("mask_cell", 1, nghost);

    // Init vof and iblanks
    init_vof_only(vof);
    init_iblank_node(iblank_node);
    init_iblank_cell(iblank_cell);

    // Create masks
    kynema_sgf::overset_ops::iblank_node_to_mask_vof(
        iblank_node, vof, mask_node);
    kynema_sgf::overset_ops::prepare_mask_cell_for_mac(repo);

    // Check against expectations
    amrex::Real error_node = check_iblank_node_impl(mask_node);
    amrex::Real error_cell = check_iblank_cell_impl(mask_cell);
    amrex::ParallelDescriptor::ReduceRealSum(error_node);
    amrex::ParallelDescriptor::ReduceRealSum(error_cell);
    EXPECT_NEAR(
        error_node, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_NEAR(
        error_cell, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // Change mask_cell to default (single-phase) approach
    kynema_sgf::overset_ops::revert_mask_cell_after_mac(repo);

    // Check against expectations
    error_cell = check_iblank_cell_default_impl(mask_cell);
    amrex::ParallelDescriptor::ReduceRealSum(error_cell);
    EXPECT_NEAR(
        error_cell, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
}

} // namespace amr_wind_tests
