#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/MassInflowOutflowBC.H"
#include "src/utilities/index_operations.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

// Count cells in box where value != offset + 1 + 20 j + 30 k (x-normal tests)
amrex::Long count_pattern_mismatches(
    const amrex::MultiFab& mfab,
    const amrex::Box& check_box,
    const amrex::Real offset)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto bx = check_box & mfab[mfi].box();
        if (!bx.ok()) {
            continue;
        }
        const auto& arr = mfab.const_array(mfi);
        reduce_op.eval(
            bx, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                const amrex::Real val =
                    offset + 1.0_rt + (20.0_rt * j) + (30.0_rt * k);
                return {static_cast<amrex::Long>(
                    std::abs(arr(i, j, k) - val) > 1.0e-12_rt)};
            });
    }
    auto nbad = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceLongSum(nbad);
    return nbad;
}

void set_pattern_and_extend(
    amrex::MultiFab& mfab,
    const amrex::Box& fill_box,
    const amrex::Orientation& ori)
{
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto bx = fill_box & mfab[mfi].box();
        if (!bx.ok()) {
            continue;
        }
        const auto& arr = mfab.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i, j, k) = 1.0_rt + (20.0_rt * j) + (30.0_rt * k);
        });
        kynema_sgf::utils::extend_boundary_ghost_layers(bx, arr, ori, 0, 1);
    }
}

void initialize_interior_pattern(amrex::MultiFab& mfab)
{
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.validbox();
        const auto& arr = mfab.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i, j, k) =
                1.0_rt + (20.0_rt * j) + (30.0_rt * k) + (1000.0_rt * i);
        });
    }
}

} // namespace

class GhostLayerFillTest : public MeshTest
{};

TEST_F(GhostLayerFillTest, extend_cell_centered_low)
{
    populate_parameters();
    initialize_mesh();

    constexpr int ng = 3;
    constexpr int idir = 0;
    const auto& domain = mesh().Geom(0).Domain();
    amrex::MultiFab mfab(mesh().boxArray(0), mesh().DistributionMap(0), 1, ng);
    mfab.setVal(0.0_rt);

    const amrex::Orientation ori(idir, amrex::Orientation::low);
    set_pattern_and_extend(mfab, amrex::adjCellLo(domain, idir, 1), ori);

    // Every ghost layer at xlo matches the first ghost layer
    EXPECT_EQ(
        count_pattern_mismatches(
            mfab, amrex::adjCellLo(domain, idir, ng), 0.0_rt),
        0);
    // Valid cells are untouched
    EXPECT_EQ(mfab.min(0, 0), 0.0_rt);
    EXPECT_EQ(mfab.max(0, 0), 0.0_rt);
}

TEST_F(GhostLayerFillTest, extend_face_normal_high)
{
    populate_parameters();
    initialize_mesh();

    constexpr int ng = 3;
    constexpr int idir = 0;
    const auto& domain = mesh().Geom(0).Domain();
    const auto ba = amrex::convert(
        mesh().boxArray(0), amrex::IntVect::TheDimensionVector(idir));
    amrex::MultiFab mfab(ba, mesh().DistributionMap(0), 1, ng);
    mfab.setVal(0.0_rt);

    // Nodal in x: boundary face at bigEnd + 1, ghost faces beyond it
    const amrex::Orientation ori(idir, amrex::Orientation::high);
    const int bnd_face = domain.bigEnd(idir) + 1;
    amrex::Box first_layer = amrex::surroundingNodes(domain, idir);
    first_layer.setRange(idir, bnd_face + 1);
    set_pattern_and_extend(mfab, first_layer, ori);

    amrex::Box all_layers = first_layer;
    all_layers.setRange(idir, bnd_face + 1, ng);
    EXPECT_EQ(count_pattern_mismatches(mfab, all_layers, 0.0_rt), 0);

    // The boundary face itself is not modified
    amrex::Box bnd_box = first_layer;
    bnd_box.setRange(idir, bnd_face);
    EXPECT_EQ(
        count_pattern_mismatches(mfab, bnd_box, 0.0_rt), bnd_box.numPts());
}

TEST_F(GhostLayerFillTest, mass_inflow_outflow_all_ghost_layers)
{
    populate_parameters();
    initialize_mesh();

    constexpr int ng = 3;
    constexpr int idir = 0;
    constexpr amrex::Real sentinel = -99.0_rt;
    auto& repo = mesh().field_repo();
    auto& velocity = repo.declare_field("velocity", 3, ng);
    auto& scalar = repo.declare_field("scalar", 1, ng);
    const auto& domain = mesh().Geom(0).Domain();

    // Flow in +x everywhere: xhi is outflow, xlo is inflow
    velocity.setVal(1.0_rt, 0, 1, ng);
    velocity.setVal(0.0_rt, 1, 2, ng);

    // Interior pattern offset by 1000 i
    scalar.setVal(sentinel);
    initialize_interior_pattern(scalar(0));

    for (const auto side :
         {amrex::Orientation::low, amrex::Orientation::high}) {
        kynema_sgf::MassInflowOutflowBC bc(
            scalar, amrex::Orientation(idir, side));
        bc(scalar, kynema_sgf::FieldState::New);
    }

    // xhi outflow: every ghost layer holds the adjacent interior value
    const amrex::Real hi_offset = 1000.0_rt * domain.bigEnd(idir);
    EXPECT_EQ(
        count_pattern_mismatches(
            scalar(0), amrex::adjCellHi(domain, idir, ng), hi_offset),
        0);

    // xlo inflow: ghost layers are left untouched
    const auto lo_layers = amrex::adjCellLo(domain, idir, ng);
    EXPECT_EQ(scalar(0).min(lo_layers, 0, ng), sentinel);
    EXPECT_EQ(scalar(0).max(lo_layers, 0, ng), sentinel);
}

} // namespace kynema_sgf_tests
