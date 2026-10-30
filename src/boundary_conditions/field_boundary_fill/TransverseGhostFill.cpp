#include "src/boundary_conditions/field_boundary_fill/TransverseGhostFill.H"
#include "src/core/FieldBCOps.H"

namespace kynema_sgf {

using namespace amrex::literals;

void fill_transverse_ghosts(
    amrex::MultiFab& mfab,
    const amrex::Geometry& geom,
    const amrex::Vector<amrex::BCRec>& bcrec,
    const amrex::Orientation ori,
    const int dcomp,
    const int ncomp,
    const int orig_comp)
{
    const int dir = ori.coordDir();
    const auto ng = mfab.nGrowVect();
    if (ng[dir] == 0) {
        return;
    }
    AMREX_ALWAYS_ASSERT(orig_comp + ncomp <= static_cast<int>(bcrec.size()));
    const auto field_domain = amrex::convert(geom.Domain(), mfab.ixType());
    auto domain = geom.Domain();
    domain.setBig(dir, field_domain.bigEnd(dir));
    auto face = domain;
    face.setRange(
        dir,
        ori.isLow() ? domain.smallEnd(dir) - ng[dir] : domain.bigEnd(dir) + 1,
        ng[dir]);

    // Auxiliary source planes are cell-centered transversely, including when
    // filling MAC fields. Complete their upper transverse nodes as ghosts too.
    amrex::BoxList boxes;
    amrex::Vector<int> owners;
    amrex::Vector<int> source_indices;
    for (int idx = 0; idx < mfab.boxArray().size(); ++idx) {
        auto box = mfab.boxArray()[idx];
        if (ori.isLow() ? box.smallEnd(dir) != field_domain.smallEnd(dir)
                        : box.bigEnd(dir) != field_domain.bigEnd(dir)) {
            continue;
        }
        box.enclosedCells();
        box.setRange(dir, face.smallEnd(dir), face.length(dir));
        if (box.ok()) {
            boxes.push_back(box);
            owners.push_back(mfab.DistributionMap()[idx]);
            source_indices.push_back(idx);
        }
    }
    if (owners.empty()) {
        return;
    }
    auto transverse_ng = ng + mfab.ixType().toIntVect();
    transverse_ng[dir] = 0;
    amrex::MultiFab slab(
        amrex::BoxArray(std::move(boxes)), amrex::DistributionMapping(owners),
        ncomp, transverse_ng);
    for (amrex::MFIter mfi(slab); mfi.isValid(); ++mfi) {
        const auto src = mfab[source_indices[mfi.index()]].const_array();
        const auto dst = slab.array(mfi);
        auto source_box = mfab[source_indices[mfi.index()]].box();
        source_box.setType(amrex::IndexType::TheCellType());
        const auto box = slab[mfi].box() & source_box;
        amrex::ParallelFor(
            box, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                dst(i, j, k, n) = src(i, j, k, dcomp + n);
            });
    }
    slab.FillBoundary(geom.periodicity());

    auto extended_domain = geom.Domain();
    extended_domain.grow(dir, ng[dir] + mfab.ixType().toIntVect()[dir]);
    auto extended_geom = geom;
    extended_geom.Domain(extended_domain);
    auto transverse_bc = amrex::Vector<amrex::BCRec>(
        bcrec.begin() + orig_comp, bcrec.begin() + orig_comp + ncomp);
    for (auto& bc : transverse_bc) {
        bc.setLo(dir, amrex::BCType::int_dir);
        bc.setHi(dir, amrex::BCType::int_dir);
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            const auto other = oit();
            if (other.coordDir() == dir) {
                continue;
            }
            const int type = other.isLow() ? bc.lo(other.coordDir())
                                           : bc.hi(other.coordDir());
            // A lower-priority prescribed face must not mask this face's data.
            if (int(other) < int(ori) &&
                (type == amrex::BCType::ext_dir ||
                 type == amrex::BCType::direction_dependent)) {
                bc.set(other, amrex::BCType::foextrap);
            }
        }
    }
    amrex::GpuBndryFuncFab<FieldBCNoOp> boundary_fill{FieldBCNoOp{}};
    for (amrex::MFIter mfi(slab); mfi.isValid(); ++mfi) {
        boundary_fill(
            slab[mfi].box(), slab[mfi], 0, ncomp, extended_geom, 0.0_rt,
            transverse_bc, 0, 0);
    }

    amrex::BoxList destination_boxes;
    amrex::Vector<int> destination_owners;
    amrex::Vector<int> destination_indices;
    for (int idx = 0; idx < mfab.boxArray().size(); ++idx) {
        auto box = amrex::grow(mfab.boxArray()[idx], ng);
        box.setType(amrex::IndexType::TheCellType());
        auto extent = box;
        extent.setRange(dir, face.smallEnd(dir), face.length(dir));
        box &= extent;
        if (box.ok()) {
            destination_boxes.push_back(box);
            destination_owners.push_back(mfab.DistributionMap()[idx]);
            destination_indices.push_back(idx);
        }
    }
    amrex::MultiFab completed(
        amrex::BoxArray(std::move(destination_boxes)),
        amrex::DistributionMapping(destination_owners), ncomp, 0);
    completed.ParallelCopy(
        slab, 0, 0, ncomp, transverse_ng, amrex::IntVect(0),
        geom.periodicity());
    for (amrex::MFIter mfi(completed); mfi.isValid(); ++mfi) {
        const auto dst = mfab[destination_indices[mfi.index()]].array();
        const auto src = completed.const_array(mfi);
        const auto box = mfi.validbox();
        for (int n = 0; n < ncomp; ++n) {
            const auto bc = bcrec[orig_comp + n];
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                for (int side = int(ori) + 1; side < 2 * AMREX_SPACEDIM;
                     ++side) {
                    const amrex::Orientation other(
                        side % AMREX_SPACEDIM, side < AMREX_SPACEDIM
                                                   ? amrex::Orientation::low
                                                   : amrex::Orientation::high);
                    const int d = other.coordDir();
                    const int type = other.isLow() ? bc.lo(d) : bc.hi(d);
                    const bool outside = other.isLow()
                                             ? iv[d] < domain.smallEnd(d)
                                             : iv[d] > domain.bigEnd(d);
                    if (outside &&
                        (type == amrex::BCType::ext_dir ||
                         type == amrex::BCType::direction_dependent)) {
                        return;
                    }
                }
                dst(iv, dcomp + n) = src(iv, n);
            });
        }
    }
}

} // namespace kynema_sgf
