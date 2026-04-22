#include "amr-wind/core/ScratchField.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_Gpu.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Geometry.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

namespace {
struct SFBCNoOp
{
    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& /* iv */,
        amrex::Array4<amrex::Real> const& /* field */,
        const int /* dcomp */,
        const int /* numcomp */,
        amrex::GeometryData const& /* geom */,
        const amrex::Real /* time */,
        const amrex::BCRec* /* bcr */,
        const int /* bcomp */,
        const int /* orig_comp */) const
    {}
};

amrex::Vector<amrex::BCRec>
scratch_field_bcrec(const amrex::Geometry& geom, const int ncomp)
{
    amrex::Vector<amrex::BCRec> bcrec(ncomp);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto bcrv = geom.isPeriodic(dir) ? amrex::BCType::int_dir
                                         : amrex::BCType::hoextrap;

        for (int n = 0; n < ncomp; ++n) {
            bcrec[n].setLo(dir, bcrv);
            bcrec[n].setHi(dir, bcrv);
        }
    }

    return bcrec;
}

} // namespace

void ScratchField::fillpatch(const amrex::Real time)
{
    fillpatch(time, num_grow());
}

void ScratchField::fillpatch(const amrex::Real time, const amrex::IntVect& ng)
{
    const int nlevels = repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        fillpatch(lev, time, this->operator()(lev), ng);
    }
}

void ScratchField::fillpatch(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost)
{
    const auto& mesh = repo().mesh();
    auto bcrec = scratch_field_bcrec(mesh.Geom(lev), num_comp());
    fillpatch(lev, time, mfab, nghost, bcrec);
}

void ScratchField::fillpatch(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& nghost,
    amrex::Vector<amrex::BCRec>& bcrec)
{
    const auto& mesh = repo().mesh();
    amrex::Interpolater* mapper = &amrex::cell_cons_interp;

    if (lev == 0) {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<SFBCNoOp>> physbc(
            mesh.Geom(lev), bcrec, SFBCNoOp{});

        amrex::FillPatchSingleLevel(
            mfab, nghost, time, {&this->operator()(lev)}, {time}, 0, 0,
            num_comp(), mesh.Geom(lev), physbc, 0);
    } else {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<SFBCNoOp>> cphysbc(
            mesh.Geom(lev - 1), bcrec, SFBCNoOp{});

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<SFBCNoOp>> fphysbc(
            mesh.Geom(lev), bcrec, SFBCNoOp{});

        amrex::FillPatchTwoLevels(
            mfab, nghost, time, {&this->operator()(lev - 1)}, {time},
            {&this->operator()(lev)}, {time}, 0, 0, num_comp(),
            mesh.Geom(lev - 1), mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
            mesh.refRatio(lev - 1), mapper, bcrec, 0);
    }
}

} // namespace amr_wind
