#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/MassInflowOutflowBC.H"
#include "src/boundary_conditions/field_boundary_fill/TransverseGhostFill.H"
#include "src/boundary_conditions/field_boundary_fill/ModulatedPowerLaw.H"
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

amrex::Long count_outflow_edge_mismatches(
    const amrex::MultiFab& mfab, const amrex::Box& domain)
{
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto arr = mfab.const_array(mfi);
        reduce_op.eval(
            mfab[mfi].box(), reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                if (i <= domain.bigEnd(0)) {
                    return {0};
                }
                const int jj = amrex::min(
                    domain.bigEnd(1), amrex::max(domain.smallEnd(1), j));
                const int kk = amrex::min(
                    domain.bigEnd(2), amrex::max(domain.smallEnd(2), k));
                const amrex::Real expected = 1.0_rt + 20.0_rt * jj +
                                             30.0_rt * kk +
                                             1000.0_rt * domain.bigEnd(0);
                return {static_cast<amrex::Long>(
                    std::abs(arr(i, j, k) - expected) > 1.0e-12_rt)};
            });
    }
    auto errors = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceLongSum(errors);
    return errors;
}

amrex::Long transverse_fill_errors(
    const amrex::Geometry& geom,
    const amrex::Orientation ori,
    const amrex::IntVect& nodal,
    const int transverse_type,
    const bool mixed = false,
    const int grid_size = 4)
{
    constexpr int ng = 3;
    constexpr amrex::Real sentinel = -999.0_rt;
    const int dir = ori.coordDir();
    auto domain = geom.Domain();
    domain.setBig(dir, domain.bigEnd(dir) + nodal[dir]);
    auto ba = amrex::BoxArray(geom.Domain());
    ba.maxSize(grid_size);
    ba.convert(nodal);
    amrex::MultiFab mfab(ba, amrex::DistributionMapping(ba), 2, ng);
    mfab.setVal(sentinel);
    const int first =
        ori.isLow() ? domain.smallEnd(dir) - 1 : domain.bigEnd(dir) + 1;
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto valid = mfi.validbox();
        if (ori.isLow() ? valid.smallEnd(dir) != domain.smallEnd(dir)
                        : valid.bigEnd(dir) != domain.bigEnd(dir)) {
            continue;
        }
        auto plane = amrex::enclosedCells(valid);
        plane.setRange(dir, first);
        const auto arr = mfab.array(mfi);
        amrex::ParallelFor(plane, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i, j, k, 1) = 10.0_rt + 2.0_rt * i + 3.0_rt * j + 5.0_rt * k;
        });
        kynema_sgf::utils::extend_boundary_ghost_layers(plane, arr, ori, 1, 1);
    }
    amrex::Vector<amrex::BCRec> bcrec(3);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bcrec[2].setLo(
            d, geom.isPeriodic(d) ? amrex::BCType::int_dir : transverse_type);
        bcrec[2].setHi(
            d, geom.isPeriodic(d) ? amrex::BCType::int_dir : transverse_type);
    }
    if (mixed) {
        const int d1 = (dir + 1) % AMREX_SPACEDIM;
        const int d2 = (dir + 2) % AMREX_SPACEDIM;
        bcrec[2].setLo(d1, amrex::BCType::hoextrap);
        bcrec[2].setHi(d1, amrex::BCType::reflect_odd);
        bcrec[2].setLo(d2, amrex::BCType::ext_dir);
        bcrec[2].setHi(d2, amrex::BCType::foextrap);
    }
    bcrec[2].set(ori, amrex::BCType::ext_dir);
    kynema_sgf::fill_transverse_ghosts(mfab, geom, bcrec, ori, 1, 1, 2);

    const auto geometry = geom.data();
    const auto bc = bcrec[2];
    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto arr = mfab.const_array(mfi);
        reduce_op.eval(
            mfab[mfi].box(), reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                const amrex::IntVect iv(i, j, k);
                const bool normal_ghost = ori.isLow()
                                              ? iv[dir] < domain.smallEnd(dir)
                                              : iv[dir] > domain.bigEnd(dir);
                amrex::Real expected = sentinel;
                if (normal_ghost) {
                    amrex::GpuArray<amrex::Real, 3> sample{
                        amrex::Real(i), amrex::Real(j), amrex::Real(k)};
                    sample[dir] = first;
                    amrex::Real sign = 1.0_rt;
                    bool preserve = false;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        if (d == dir) {
                            continue;
                        }
                        const int lo = domain.smallEnd(d);
                        const int hi = domain.bigEnd(d);
                        if (geometry.isPeriodic(d)) {
                            const int length = geometry.Domain().length(d);
                            if (iv[d] < lo) {
                                sample[d] += length;
                            } else if (iv[d] > hi) {
                                sample[d] -= length;
                            }
                            continue;
                        }
                        const bool low = iv[d] < lo;
                        const bool high = iv[d] > hi;
                        if (!low && !high) {
                            continue;
                        }
                        const amrex::Orientation other(
                            d, low ? amrex::Orientation::low
                                   : amrex::Orientation::high);
                        const int type = low ? bc.lo(d) : bc.hi(d);
                        if (type == amrex::BCType::ext_dir ||
                            type == amrex::BCType::direction_dependent) {
                            preserve = preserve || int(other) > int(ori);
                            sample[d] = low ? lo : hi;
                        } else if (
                            type == amrex::BCType::reflect_even ||
                            type == amrex::BCType::reflect_odd) {
                            sample[d] =
                                low ? 2 * lo - iv[d] - 1 : 2 * hi - iv[d] + 1;
                            if (type == amrex::BCType::reflect_odd) {
                                sign = -sign;
                            }
                        } else {
                            sample[d] = low ? lo : hi;
                            if (type == amrex::BCType::hoextrap &&
                                (low ? iv[d] == lo - 1 : iv[d] == hi + 1)) {
                                sample[d] += low ? -0.5_rt : 0.5_rt;
                            }
                        }
                    }
                    if (!preserve) {
                        expected =
                            sign * (10.0_rt + 2.0_rt * sample[0] +
                                    3.0_rt * sample[1] + 5.0_rt * sample[2]);
                    }
                }
                return {static_cast<amrex::Long>(
                    std::abs(arr(iv, 1) - expected) > 1.0e-12_rt ||
                    arr(iv, 0) != sentinel)};
            });
    }
    auto errors = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceLongSum(errors);
    return errors;
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
    amrex::ParmParse("geometry")
        .addarr("is_periodic", amrex::Vector<int>{0, 0, 0});
    initialize_mesh();

    constexpr int ng = 3;
    constexpr int idir = 0;
    constexpr amrex::Real sentinel = -99.0_rt;
    auto& repo = mesh().field_repo();
    auto& velocity = repo.declare_field("velocity", 3, ng);
    auto& scalar = repo.declare_field("scalar", 1, ng);
    const auto& domain = mesh().Geom(0).Domain();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        scalar.bcrec()[0].setLo(d, amrex::BCType::foextrap);
        scalar.bcrec()[0].setHi(d, amrex::BCType::foextrap);
    }

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
    EXPECT_EQ(count_outflow_edge_mismatches(scalar(0), domain), 0);

    // xlo inflow: ghost layers are left untouched
    const auto lo_layers = amrex::adjCellLo(domain, idir, ng);
    EXPECT_EQ(scalar(0).min(lo_layers, 0, ng), sentinel);
    EXPECT_EQ(scalar(0).max(lo_layers, 0, ng), sentinel);
}

TEST_F(GhostLayerFillTest, transverse_physical_and_grid_boundaries)
{
    const amrex::Box domain(
        amrex::IntVect(4, 6, 8), amrex::IntVect(11, 13, 15));
    const amrex::RealBox physical_box(
        {0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt});
    const amrex::Geometry geom(
        domain, physical_box, 0, amrex::Array<int, 3>{0, 0, 0});
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        EXPECT_EQ(
            transverse_fill_errors(
                geom, ori, amrex::IntVect(0), amrex::BCType::foextrap, true),
            0);
        EXPECT_EQ(
            transverse_fill_errors(
                geom, ori, amrex::IntVect(0), amrex::BCType::foextrap, false,
                2),
            0);
        for (const int type :
             {amrex::BCType::foextrap, amrex::BCType::hoextrap,
              amrex::BCType::reflect_even, amrex::BCType::reflect_odd,
              amrex::BCType::ext_dir, amrex::BCType::direction_dependent}) {
            SCOPED_TRACE(
                ::testing::Message() << "face " << int(ori) << ", BC " << type);
            EXPECT_EQ(
                transverse_fill_errors(geom, ori, amrex::IntVect(0), type), 0);
            EXPECT_EQ(
                transverse_fill_errors(
                    geom, ori,
                    amrex::IntVect::TheDimensionVector(ori.coordDir()), type),
                0);
            EXPECT_EQ(
                transverse_fill_errors(
                    geom, ori,
                    amrex::IntVect::TheDimensionVector(
                        (ori.coordDir() + 1) % AMREX_SPACEDIM),
                    type),
                0);
        }
    }
}

TEST_F(GhostLayerFillTest, transverse_periodic_boundaries)
{
    const amrex::Box domain(
        amrex::IntVect(4, 6, 8), amrex::IntVect(11, 13, 15));
    const amrex::RealBox physical_box(
        {0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt});
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        for (int mask = 1; mask < 4; ++mask) {
            amrex::Array<int, 3> periodic{0, 0, 0};
            periodic[(ori.coordDir() + 1) % AMREX_SPACEDIM] = mask & 1;
            periodic[(ori.coordDir() + 2) % AMREX_SPACEDIM] = (mask >> 1) & 1;
            const amrex::Geometry geom(domain, physical_box, 0, periodic);
            EXPECT_EQ(
                transverse_fill_errors(
                    geom, ori, amrex::IntVect(0), amrex::BCType::foextrap),
                0);
            EXPECT_EQ(
                transverse_fill_errors(
                    geom, ori,
                    amrex::IntVect::TheDimensionVector(
                        (ori.coordDir() + 1) % AMREX_SPACEDIM),
                    amrex::BCType::foextrap),
                0);
        }
    }
}

TEST_F(GhostLayerFillTest, modulated_power_law_completes_transverse_ghosts)
{
    populate_parameters();
    amrex::ParmParse("geometry")
        .addarr("is_periodic", amrex::Vector<int>{0, 0, 0});
    amrex::ParmParse("ABL").addarr(
        "temperature_heights", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt});
    amrex::ParmParse("ABL").addarr(
        "temperature_values", amrex::Vector<amrex::Real>{300.0_rt, 300.0_rt});
    amrex::ParmParse("MPL").add("delta_t", 0.0_rt);
    initialize_mesh();
    auto& repo = sim().repo();
    auto& velocity = repo.declare_field("velocity", 3, 3);
    auto& temperature = repo.declare_field("temperature", 1, 3);
    const amrex::Orientation ori(0, amrex::Orientation::low);
    for (auto* field : {&velocity, &temperature}) {
        field->setVal(-99.0_rt);
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            field->bc_type()[oit()] = BC::zero_gradient;
        }
        field->bc_type()[ori] = BC::mass_inflow;
        for (auto& bc : field->bcrec()) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                bc.setLo(d, amrex::BCType::foextrap);
                bc.setHi(d, amrex::BCType::foextrap);
            }
            bc.set(ori, amrex::BCType::ext_dir);
        }
    }
    kynema_sgf::ModulatedPowerLaw profile(sim());
    profile.set_temperature(0, 0.0_rt, temperature, temperature(0));
    profile.set_velocity(0, 0.0_rt, velocity, velocity(0));
    auto face = amrex::grow(mesh().Geom(0).Domain(), 3);
    face.setRange(0, mesh().Geom(0).Domain().smallEnd(0) - 3, 3);
    EXPECT_EQ(temperature(0).min(face, 0, 3), 300.0_rt);
    EXPECT_EQ(temperature(0).max(face, 0, 3), 300.0_rt);
    for (int n = 0; n < 3; ++n) {
        EXPECT_GE(velocity(0).min(face, n, 3), 0.0_rt);
    }
    EXPECT_EQ(temperature(0).min(0, 0), -99.0_rt);
    EXPECT_EQ(temperature(0).max(0, 0), -99.0_rt);
}

} // namespace kynema_sgf_tests
