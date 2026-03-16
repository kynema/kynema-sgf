#ifndef FIELDFILLPATCHOPS_H
#define FIELDFILLPATCHOPS_H

#include <utility>

#include "amr-wind/core/Field.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldBCOps.H"

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_REAL.H"
#include "AMReX_PhysBCFunct.H"
#include "AMReX_FillPatchUtil.H"

/**
 *  \defgroup field_fillpatch Fillpatch operators
 *  Field fillpatch operators
 *
 *  This group documents classes and functions available for performing
 *  fillpatch operations on the entire field. This API represents a higher-level
 *  interface to the MultiFab based fillpatch operations and fill physics BC
 *  operations defined in AMReX library. The functionality available in AMR-Wind
 *  has two components: 1. the fill patch operators themselves that perform the
 *  fill operation, and 2. operator classes that define what happens during a
 *  fillpatch operation at the domain boundaries.
 *
 *  \ingroup fields
 */

namespace amr_wind {

/** Abstract interface for classes that handle fill patch and physical boundary
 *  condition operations for multi-level AMR meshes.
 *
 *  \ingroup field_fillpatch
 */
class FieldFillPatchOpsBase
{
public:
    FieldFillPatchOpsBase() = default;

    virtual ~FieldFillPatchOpsBase() = default;

    //! Implementation that handles filling patches on a single level as well as
    //! across a coarse-fine interface.
    virtual void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) = 0;

    //! Implementation that handles filling patches on a single level as well as
    //! across a coarse-fine interface, for an array of different fields
    virtual void fillpatch_sibling_fields(
        int lev,
        amrex::Real time,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
        const amrex::IntVect& nghost,
        const amrex::Vector<amrex::BCRec>& fillpatch_bcrec,
        const amrex::Vector<amrex::BCRec>& physbc_bcrec,
        FieldState fstate = FieldState::New) = 0;

    //! Implementation that handles filling patches from a coarse to fine level
    virtual void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) = 0;

    //! Implementation that handles filling physical boundary conditions
    virtual void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) = 0;

    virtual void set_inflow(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) = 0;

    virtual void set_inflow_sibling_fields(
        int lev,
        amrex::Real time,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs) = 0;
};

/** Implementation that just fills a constant value on newly created grids
 *
 *  \ingroup field_fillpatch
 */
class FieldFillConstScalar : public FieldFillPatchOpsBase
{
public:
    FieldFillConstScalar(Field& /*unused*/, amrex::Real fill_val)
        : m_fill_val(fill_val)
    {}

    void fillpatch(
        int /*lev*/,
        amrex::Real /*time*/,
        amrex::MultiFab& mfab,
        const amrex::IntVect& /*nghost*/,
        const FieldState /*fstate*/) override
    {
        mfab.setVal(m_fill_val);
    }

    void fillpatch_sibling_fields(
        int /*lev*/,
        amrex::Real /*time*/,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& /*ffabs*/,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& /*cfabs*/,
        const amrex::IntVect& /*nghost*/,
        const amrex::Vector<amrex::BCRec>& /*bcrec*/,
        const amrex::Vector<amrex::BCRec>& /*bcrec*/,
        const FieldState /*fstate*/) override
    {
        for (const auto& mfab : mfabs) {
            mfab->setVal(m_fill_val);
        }
    }

    void fillpatch_from_coarse(
        int /*lev*/,
        amrex::Real /*time*/,
        amrex::MultiFab& mfab,
        const amrex::IntVect& /*nghost*/,
        const FieldState /*fstate*/) override
    {
        mfab.setVal(m_fill_val);
    }

    void fillphysbc(
        int /*lev*/,
        amrex::Real /*time*/,
        amrex::MultiFab& mfab,
        const amrex::IntVect& /*nghost*/,
        const FieldState /*fstate*/) override
    {
        mfab.setVal(m_fill_val);
    }

    void set_inflow(
        int /*lev*/,
        amrex::Real /*time*/,
        amrex::MultiFab& /*mfab*/,
        const amrex::IntVect& /*nghost*/,
        const FieldState /*fstate*/) override
    {
        amrex::Abort("FieldFillConstScalar::set_inflow is not implemented");
    }

private:
    amrex::Real m_fill_val;
};

/** Concrete implementation of filling patches and boundary conditions on
 *  multi-level AMR meshes.
 *
 *  \ingroup field_fillpatch
 */
template <typename BCOpCreator>
class FieldFillPatchOps : public FieldFillPatchOpsBase
{
public:
    using Functor = typename BCOpCreator::FunctorType;
    /**
     *  \param field Field whose patches are filled by this instance
     *  \param mesh  The mesh instance to determine amrex::Geometry at a level
     *  \param time  Time instance
     *  \param itype Interpolation type (default: conservative cell linear)
     *  \param face_itype Face interpolation type (default: face divergence
     * free)
     */
    FieldFillPatchOps(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        FieldInterpolator itype = FieldInterpolator::CellConsLinear,
        FieldInterpolator face_itype = FieldInterpolator::FaceDivFree)
        : m_time(time)
        , m_mesh(mesh)
        , m_field(field)
        , m_op(field)
        , m_mapper(field_impl::get_interpolation_operator(itype))
        , m_face_mapper(field_impl::get_interpolation_operator(face_itype))
    {
        check_face_mapper();
    }

    FieldFillPatchOps(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        const BCOpCreator& bc_op,
        FieldInterpolator itype = FieldInterpolator::CellConsLinear,
        FieldInterpolator face_itype = FieldInterpolator::FaceDivFree)
        : m_time(time)
        , m_mesh(mesh)
        , m_field(field)
        , m_op(bc_op)
        , m_mapper(field_impl::get_interpolation_operator(itype))
        , m_face_mapper(field_impl::get_interpolation_operator(face_itype))
    {
        check_face_mapper();
    }

    /** Utility function to get a vector of multifabs for over states
     *
     *  Returns one or two multifabs at a given level depending on whether the
     *  field was registered with states. When two states are returned the first
     *  one is the Old state and the second one is the New state.
     */
    amrex::Vector<amrex::MultiFab*> get_mfab_vec(int lev)
    {
        const int nstates = amrex::min(m_field.num_time_states(), 2);
        amrex::Vector<amrex::MultiFab*> ret;

        // The states in the FieldInfo data are ordered from newest to oldest,
        // so swap the order
        for (int i = nstates - 1; i >= 0; --i) {
            const auto fstate = static_cast<FieldState>(i);
            ret.push_back(&m_field.state(fstate)(lev));
        }
        return ret;
    }

    // Version that does no interpolation in time

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState fstate = FieldState::New) override
    {
        auto& fld = m_field.state(fstate);
        if (lev == 0) {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
                m_mesh.Geom(lev), m_field.bcrec(), bc_functor());

            amrex::FillPatchSingleLevel(
                mfab, nghost, time, {&fld(lev)}, {time}, 0, 0,
                m_field.num_comp(), m_mesh.Geom(lev), physbc, 0);
        } else {
            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
                m_mesh.Geom(lev - 1), m_field.bcrec(), bc_functor());

            amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
                m_mesh.Geom(lev), m_field.bcrec(), bc_functor());

            amrex::FillPatchTwoLevels(
                mfab, nghost, time, {&fld(lev - 1)}, {time}, {&fld(lev)},
                {time}, 0, 0, m_field.num_comp(), m_mesh.Geom(lev - 1),
                m_mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
                m_mesh.refRatio(lev - 1), m_mapper, m_field.bcrec(), 0);
        }
    }

    void fillpatch_sibling_fields(
        int lev,
        amrex::Real time,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
        const amrex::IntVect& nghost,
        const amrex::Vector<amrex::BCRec>& fillpatch_bcrec,
        const amrex::Vector<amrex::BCRec>& physbc_bcrec,
        const FieldState /*fstate = FieldState::New*/) override
    {
        if (lev == 0) {
            for (int i = 0; std::cmp_less(i, mfabs.size()); i++) {
                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
                    m_mesh.Geom(lev), physbc_bcrec, bc_functor_face(i));
                amrex::FillPatchSingleLevel(
                    *mfabs[i], nghost, time, {ffabs[i]}, {time}, 0, 0, 1,
                    m_mesh.Geom(lev), physbc, i);
            }
        } else {
            AMREX_D_TERM(
                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbcx(
                    m_mesh.Geom(lev - 1), physbc_bcrec, bc_functor_face(0));
                , amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbcy(
                      m_mesh.Geom(lev - 1), physbc_bcrec, bc_functor_face(1));
                , amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbcz(
                      m_mesh.Geom(lev - 1), physbc_bcrec, bc_functor_face(2)););

            AMREX_D_TERM(
                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbcx(
                    m_mesh.Geom(lev), physbc_bcrec, bc_functor_face(0));
                , amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbcy(
                      m_mesh.Geom(lev), physbc_bcrec, bc_functor_face(1));
                , amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbcz(
                      m_mesh.Geom(lev), physbc_bcrec, bc_functor_face(2)););

            amrex::Array<
                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>>,
                AMREX_SPACEDIM>
                cphysbc_arr = {AMREX_D_DECL(cphysbcx, cphysbcy, cphysbcz)};

            amrex::Array<
                amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>>,
                AMREX_SPACEDIM>
                fphysbc_arr = {AMREX_D_DECL(fphysbcx, fphysbcy, fphysbcz)};

            amrex::Array<int, AMREX_SPACEDIM> idx = {AMREX_D_DECL(0, 1, 2)};
            const amrex::Array<amrex::Vector<amrex::BCRec>, AMREX_SPACEDIM>
                bcrec_arr = {AMREX_D_DECL(
                    fillpatch_bcrec, fillpatch_bcrec, fillpatch_bcrec)};

            amrex::FillPatchTwoLevels(
                mfabs, nghost, time, {cfabs}, {time}, {ffabs}, {time}, 0, 0, 1,
                m_mesh.Geom(lev - 1), m_mesh.Geom(lev), cphysbc_arr, idx,
                fphysbc_arr, idx, m_mesh.refRatio(lev - 1), m_face_mapper,
                bcrec_arr, idx);
        }
    }

    void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState fstate = FieldState::New) override
    {
        const auto& fld = m_field.state(fstate);
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> cphysbc(
            m_mesh.Geom(lev - 1), m_field.bcrec(), bc_functor());

        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> fphysbc(
            m_mesh.Geom(lev), m_field.bcrec(), bc_functor());

        amrex::InterpFromCoarseLevel(
            mfab, nghost, time, fld(lev - 1), 0, 0, m_field.num_comp(),
            m_mesh.Geom(lev - 1), m_mesh.Geom(lev), cphysbc, 0, fphysbc, 0,
            m_mesh.refRatio(lev - 1), m_mapper, m_field.bcrec(), 0);
    }

    void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState /*fstate*/) override
    {
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<Functor>> physbc(
            m_mesh.Geom(lev), m_field.bcrec(), bc_functor());
        physbc.FillBoundary(mfab, 0, m_field.num_comp(), nghost, time, 0);
    }

    void set_inflow(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        const FieldState /*fstate*/) override
    {
        const int ng = nghost[0];
        const auto& bctype = m_field.bc_type();
        const auto& geom = m_mesh.Geom(lev);
        const auto& gdata = geom.data();
        const auto& domain = geom.growPeriodicDomain(ng);
        const auto& bcfunc = bc_functor();
        const auto& ncomp = m_field.num_comp();

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            const auto ori = oit();
            if ((bctype[ori] != BC::mass_inflow) &&
                (bctype[ori] != BC::mass_inflow_outflow)) {
                continue;
            }

            const int idir = ori.coordDir();
            const auto& dbx =
                ori.isLow() ? amrex::adjCellLo(domain, idir, nghost[idir])
                            : amrex::adjCellHi(domain, idir, nghost[idir]);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
                const auto& gbx = amrex::grow(mfi.validbox(), nghost);
                const auto& bx = gbx & dbx;
                if (!bx.ok()) {
                    continue;
                }

                const auto& marr = mfab[mfi].array();
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(
                            const int i, const int j, const int k) {
                        for (int n = 0; n < ncomp; ++n) {
                            bcfunc.set_inflow(
                                {i, j, k}, marr, gdata, time, ori, n, 0, 0);
                        }
                    });
            }
        }
    }

    void set_inflow_sibling_fields(
        const int lev,
        const amrex::Real time,
        const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs) override
    {
        const auto& bctype = m_field.bc_type();
        const auto& geom = m_mesh.Geom(lev);
        const auto& gdata = geom.data();
        const auto& bcfunc = bc_functor();

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            const auto ori = oit();
            if ((bctype[ori] != BC::mass_inflow) &&
                (bctype[ori] != BC::mass_inflow_outflow) &&
                (bctype[ori] != BC::wave_generation)) {
                continue;
            }

            const int idir = ori.coordDir();
            const auto& domain_box = geom.Domain();
            for (int fdir = 0; fdir < AMREX_SPACEDIM; ++fdir) {

                // Only face-normal velocities populated here
                if (idir != fdir) {
                    continue;
                }
                const auto& dbx = ori.isLow() ? amrex::bdryLo(domain_box, idir)
                                              : amrex::bdryHi(domain_box, idir);

                auto& mfab = *mfabs[fdir];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.validbox();
                    const auto& bx = vbx & dbx;
                    if (!bx.ok()) {
                        continue;
                    }

                    const auto& marr = mfab[mfi].array();
                    amrex::FArrayBox tmp_fab(
                        bx, mfab.nComp(), amrex::The_Async_Arena());
                    amrex::Array4<amrex::Real> const& tmp_marr =
                        tmp_fab.array();

                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(
                                const int i, const int j, const int k) {
                            bcfunc.set_inflow(
                                {i, j, k}, tmp_marr, gdata, time, ori, 0, 0,
                                fdir);

                            const bool is_outflow =
                                (ori.isLow() && tmp_marr(i, j, k) < 0) ||
                                (ori.isHigh() && tmp_marr(i, j, k) > 0);
                            if ((bctype[ori] == BC::mass_inflow_outflow) &&
                                is_outflow) {
                                marr(i, j, k) = marr(i, j, k);
                            } else {
                                marr(i, j, k) = tmp_marr(i, j, k);
                            }
                        });
                }
            }
        }
    }

protected:
    Functor bc_functor() { return m_op(); }

    Functor bc_functor_face(const int face_dir) { return m_op(face_dir); }

    void check_face_mapper()
    {
        if (std::any_of(
                m_mesh.refRatio().begin(), m_mesh.refRatio().end(),
                [](amrex::IntVect iv) { return iv != amrex::IntVect(2); })) {
            amrex::Print()
                << "WARNING: Changing the face interpolator to FaceLinear "
                   "because a refinement ratio is different than 2"
                << '\n';
            m_face_mapper = field_impl::get_interpolation_operator(
                FieldInterpolator::FaceLinear);
        }
    }

    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    Field& m_field;

    const BCOpCreator m_op;

    //! Function that handles interpolation from coarse to fine level
    amrex::Interpolater* m_mapper;

    //! Function that handles interpolation from coarse to fine level for faces
    amrex::Interpolater* m_face_mapper;
};

} // namespace amr_wind

#endif /* FIELDFILLPATCHOPS_H */
