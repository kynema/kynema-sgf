#include "src/utilities/index_operations.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_Interpolater.H"

namespace kynema_sgf::utils {

using namespace amrex::literals;

amrex::FArrayBox interpolate_boundary_plane(
    const amrex::FArrayBox& coarse,
    const amrex::FArrayBox& fine,
    const amrex::Orientation ori,
    amrex::IntVect ratio,
    const amrex::Geometry& coarse_geom,
    const amrex::Geometry& fine_geom,
    const int nghost)
{
    const int dir = ori.coordDir();
    AMREX_ALWAYS_ASSERT(coarse.nComp() == fine.nComp());
    AMREX_ALWAYS_ASSERT(coarse.box().length(dir) == 1);
    AMREX_ALWAYS_ASSERT(fine.box().length(dir) == 1);
    auto region = fine.box();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != dir) {
            region.grow(d, nghost);
            if (!fine_geom.isPeriodic(d)) {
                region.setSmall(
                    d, amrex::max(
                           region.smallEnd(d), fine_geom.Domain().smallEnd(d)));
                region.setBig(
                    d,
                    amrex::min(region.bigEnd(d), fine_geom.Domain().bigEnd(d)));
            }
        }
    }
    amrex::FArrayBox completed(region, fine.nComp());
    if (region == fine.box()) {
        completed.copy<amrex::RunOn::Device>(fine);
        return completed;
    }

    // Interpolate only transversely: the plane normal is a single cell at
    // index zero on both levels, including for high-side boundary planes.
    ratio[dir] = 1;
    auto flat_region = region;
    flat_region.setRange(dir, 0, 1);
    const auto stencil = amrex::cell_cons_interp.CoarseBox(flat_region, ratio);
    auto domain = coarse_geom.Domain();
    domain.setRange(dir, 0, 1);
    auto required = stencil & domain;
    required.setRange(dir, coarse.box().smallEnd(dir), 1);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        coarse.box().contains(required),
        "BoundaryPlane: coarse plane does not cover the fine interpolation "
        "stencil");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != dir && coarse_geom.isPeriodic(d) &&
            (stencil.smallEnd(d) < domain.smallEnd(d) ||
             stencil.bigEnd(d) > domain.bigEnd(d))) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                coarse.box().smallEnd(d) <= domain.smallEnd(d) &&
                    coarse.box().bigEnd(d) >= domain.bigEnd(d),
                "BoundaryPlane: coarse plane does not cover a periodic "
                "interpolation stencil");
        }
    }
    amrex::FArrayBox source(stencil, coarse.nComp());
    const auto src = source.array();
    const auto data = coarse.const_array();
    const auto lo = domain.smallEnd();
    const auto hi = domain.bigEnd();
    const auto periodic = coarse_geom.isPeriodicArray();
    const int normal = coarse.box().smallEnd(dir);
    amrex::ParallelFor(
        stencil, coarse.nComp(),
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
            amrex::IntVect sample(i, j, k);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (d == dir) {
                    sample[d] = normal;
                } else if (periodic[d] == 1) {
                    const int length = hi[d] - lo[d] + 1;
                    sample[d] =
                        lo[d] +
                        ((((sample[d] - lo[d]) % length) + length) % length);
                } else {
                    sample[d] = amrex::min(hi[d], amrex::max(lo[d], sample[d]));
                }
            }
            src(i, j, k, n) = data(sample, n);
        });
    auto flat_coarse_geom = coarse_geom;
    flat_coarse_geom.Domain(domain);
    auto flat_fine_geom = fine_geom;
    auto fine_domain = fine_geom.Domain();
    fine_domain.setRange(dir, 0, 1);
    flat_fine_geom.Domain(fine_domain);
    amrex::Vector<amrex::BCRec> bcrec(fine.nComp());
    for (auto& bc : bcrec) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const int type = periodic[d] == 1 || d == dir
                                 ? amrex::BCType::int_dir
                                 : amrex::BCType::foextrap;
            bc.setLo(d, type);
            bc.setHi(d, type);
        }
    }
    amrex::FArrayBox interpolated(flat_region, fine.nComp());
    amrex::cell_cons_interp.interp(
        source, 0, interpolated, 0, fine.nComp(), flat_region, ratio,
        flat_coarse_geom, flat_fine_geom, bcrec, 0, 0, amrex::RunOn::Device);
    const auto out = completed.array();
    const auto values = interpolated.const_array();
    amrex::ParallelFor(
        region, fine.nComp(), [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
            amrex::IntVect sample(i, j, k);
            sample[dir] = 0;
            out(i, j, k, n) = values(sample, n);
        });
    const auto recorded = fine.const_array();
    const auto recorded_box = fine.box();
    const auto fine_lo = fine_domain.smallEnd();
    const auto fine_hi = fine_domain.bigEnd();
    amrex::ParallelFor(
        region, fine.nComp(), [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
            amrex::IntVect sample(i, j, k);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (d != dir && periodic[d] == 1) {
                    const int length = fine_hi[d] - fine_lo[d] + 1;
                    sample[d] =
                        fine_lo[d] +
                        ((((sample[d] - fine_lo[d]) % length) + length) %
                         length);
                }
            }
            if (recorded_box.contains(sample)) {
                out(i, j, k, n) = recorded(sample, n);
            }
        });
    amrex::Gpu::streamSynchronize();
    return completed;
}

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
