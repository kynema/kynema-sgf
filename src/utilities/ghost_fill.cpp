#include "src/utilities/index_operations.H"
#include "AMReX_PhysBCFunct.H"

namespace kynema_sgf::utils {

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
    auto transverse_ng = ng + mfab.ixType().toIntVect();
    transverse_ng[dir] = 0;
    amrex::BoxList boxes;
    amrex::BoxList donor_boxes;
    amrex::Vector<int> owners;
    amrex::Vector<int> source_indices;
    amrex::BoxList destination_boxes;
    amrex::Vector<int> destination_owners;
    amrex::Vector<int> destination_indices;
    for (int idx = 0; idx < mfab.boxArray().size(); ++idx) {
        auto box = mfab.boxArray()[idx];
        auto destination = amrex::grow(box, ng);
        destination.setType(amrex::IndexType::TheCellType());
        destination.setSmall(
            dir, amrex::max(destination.smallEnd(dir), face.smallEnd(dir)));
        destination.setBig(
            dir, amrex::min(destination.bigEnd(dir), face.bigEnd(dir)));
        if (destination.ok()) {
            destination_boxes.push_back(destination);
            destination_owners.push_back(mfab.DistributionMap()[idx]);
            destination_indices.push_back(idx);
        }
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
            // Only valid face data and its physical extensions may be donors.
            // Internal ghosts can contain stale coarse/fine fallback values.
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (d == dir || geom.isPeriodic(d)) {
                    continue;
                }
                if (box.smallEnd(d) == domain.smallEnd(d)) {
                    box.growLo(d, ng[d]);
                }
                if (box.bigEnd(d) == domain.bigEnd(d)) {
                    box.growHi(d, transverse_ng[d]);
                }
            }
            donor_boxes.push_back(box);
        }
    }
    if (owners.empty()) {
        return;
    }
    amrex::MultiFab slab(
        amrex::BoxArray(std::move(boxes)), amrex::DistributionMapping(owners),
        ncomp, transverse_ng);
    // Non-owning cell-centered views keep FAB copy box types consistent for MAC
    // fields without changing their indices or allocating additional storage.
    for (amrex::MFIter mfi(slab); mfi.isValid(); ++mfi) {
        const amrex::FArrayBox src(
            mfab.const_array(source_indices[mfi.index()]));
        const auto box = slab[mfi].box() & src.box();
        slab[mfi].copy<amrex::RunOn::Device>(src, box, dcomp, box, 0, ncomp);
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
            const int type = bc.data()[int(other)];
            // A lower-priority prescribed face must not mask this face's data.
            if (int(other) < int(ori) &&
                (type == amrex::BCType::ext_dir ||
                 type == amrex::BCType::direction_dependent)) {
                bc.set(other, amrex::BCType::foextrap);
            }
        }
    }
    amrex::GpuBndryFuncFab<amrex::FabFillNoOp> boundary_fill;
    for (amrex::MFIter mfi(slab); mfi.isValid(); ++mfi) {
        auto box = slab[mfi].box();
        // Extra transverse MAC nodes exist only on the high side. The unused
        // low-side scratch cells have no field data and must not be
        // extrapolated.
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d != dir) {
                box.growLo(d, -mfab.ixType().toIntVect()[d]);
            }
        }
        boundary_fill(
            box, slab[mfi], 0, ncomp, extended_geom, 0.0_rt, transverse_bc, 0,
            0);
    }

    amrex::MultiFab donors(
        amrex::BoxArray(std::move(donor_boxes)), slab.DistributionMap(), ncomp,
        0);
    for (amrex::MFIter mfi(donors); mfi.isValid(); ++mfi) {
        const auto& box = mfi.validbox();
        donors[mfi].copy<amrex::RunOn::Device>(
            slab[mfi], box, 0, box, 0, ncomp);
    }
    amrex::MultiFab completed(
        amrex::BoxArray(std::move(destination_boxes)),
        amrex::DistributionMapping(destination_owners), ncomp, 0);
    for (amrex::MFIter mfi(completed); mfi.isValid(); ++mfi) {
        const auto& box = mfi.validbox();
        const amrex::FArrayBox src(
            mfab.const_array(destination_indices[mfi.index()]));
        completed[mfi].copy<amrex::RunOn::Device>(
            src, box, dcomp, box, 0, ncomp);
    }
    completed.ParallelCopy(
        donors, 0, 0, ncomp, amrex::IntVect(0), amrex::IntVect(0),
        geom.periodicity());
    for (amrex::MFIter mfi(completed); mfi.isValid(); ++mfi) {
        amrex::FArrayBox dst(mfab.array(destination_indices[mfi.index()]));
        for (int n = 0; n < ncomp; ++n) {
            const auto bc = bcrec[orig_comp + n];
            auto box = mfi.validbox();
            for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
                const auto other = oit();
                const int d = other.coordDir();
                const int type = bc.data()[int(other)];
                if (int(other) <= int(ori) || d == dir ||
                    (type != amrex::BCType::ext_dir &&
                     type != amrex::BCType::direction_dependent)) {
                    continue;
                }
                if (other.isLow()) {
                    box.setSmall(
                        d, amrex::max(box.smallEnd(d), domain.smallEnd(d)));
                } else {
                    box.setBig(d, amrex::min(box.bigEnd(d), domain.bigEnd(d)));
                }
            }
            if (box.ok()) {
                dst.copy<amrex::RunOn::Device>(
                    completed[mfi], box, n, box, dcomp + n, 1);
            }
        }
    }
}

} // namespace kynema_sgf::utils
