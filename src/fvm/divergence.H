#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_Array4.H"
#include "AMReX_Geometry.H"
#include "AMReX_MFIter.H"

namespace kynema_sgf::fvm {

/** Divergence operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Divergence
{
    Divergence(FTypeOut& divphi, const FTypeIn& phi)
        : m_divphi(divphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(
            AMREX_SPACEDIM * divphi.num_comp() == phi.num_comp());
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const int ncomp = m_divphi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& divphi_arr = m_divphi(lev).array(mfi);
        const auto& phi_arr = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                amrex::Real cp1 = Stencil::c00;
                amrex::Real c = Stencil::c01;
                amrex::Real cm1 = Stencil::c02;
                divphi_arr(i, j, k, icomp) =
                    (cp1 * phi_arr(i + 1, j, k, (icomp * AMREX_SPACEDIM) + 0) +
                     c * phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 0) +
                     cm1 * phi_arr(i - 1, j, k, (icomp * AMREX_SPACEDIM) + 0)) *
                    idx[0];

                cp1 = Stencil::c10;
                c = Stencil::c11;
                cm1 = Stencil::c12;
                divphi_arr(i, j, k, icomp) +=
                    (cp1 * phi_arr(i, j + 1, k, (icomp * AMREX_SPACEDIM) + 1) +
                     c * phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 1) +
                     cm1 * phi_arr(i, j - 1, k, (icomp * AMREX_SPACEDIM) + 1)) *
                    idx[1];

                cp1 = Stencil::c20;
                c = Stencil::c21;
                cm1 = Stencil::c22;
                divphi_arr(i, j, k, icomp) +=
                    (cp1 * phi_arr(i, j, k + 1, (icomp * AMREX_SPACEDIM) + 2) +
                     c * phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 2) +
                     cm1 * phi_arr(i, j, k - 1, (icomp * AMREX_SPACEDIM) + 2)) *
                    idx[2];
            }
        });
    }

    void applyFused(const int lev) const
    {
        const int ncomp = m_divphi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();
        const amrex::GpuArray<int, AMREX_SPACEDIM> periodic{
            geom.isPeriodic(0) ? 1 : 0,
            geom.isPeriodic(1) ? 1 : 0,
            geom.isPeriodic(2) ? 1 : 0};

        auto& out_mf = m_divphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const auto sc = stencil::pick_stencil(
                    i, j, k, dlo, dhi, periodic);
                auto const& phi_arr = in_arrs[nbx];
                auto const& divphi_arr = out_arrs[nbx];
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    divphi_arr(i, j, k, icomp) =
                        (sc.c00 *
                             phi_arr(i + 1, j, k, (icomp * AMREX_SPACEDIM) + 0) +
                         sc.c01 *
                             phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 0) +
                         sc.c02 *
                             phi_arr(i - 1, j, k, (icomp * AMREX_SPACEDIM) + 0)) *
                        idx[0];

                    divphi_arr(i, j, k, icomp) +=
                        (sc.c10 *
                             phi_arr(i, j + 1, k, (icomp * AMREX_SPACEDIM) + 1) +
                         sc.c11 *
                             phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 1) +
                         sc.c12 *
                             phi_arr(i, j - 1, k, (icomp * AMREX_SPACEDIM) + 1)) *
                        idx[1];

                    divphi_arr(i, j, k, icomp) +=
                        (sc.c20 *
                             phi_arr(i, j, k + 1, (icomp * AMREX_SPACEDIM) + 2) +
                         sc.c21 *
                             phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 2) +
                         sc.c22 *
                             phi_arr(i, j, k - 1, (icomp * AMREX_SPACEDIM) + 2)) *
                        idx[2];
                }
            });
        amrex::Gpu::streamSynchronize();
    }

    FTypeOut& m_divphi;
    const FTypeIn& m_phi;
};

/** Compute the divergence of a vector field
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
void divergence(FTypeOut& divphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::divergence");
    Divergence<FTypeIn, FTypeOut> grad(divphi, phi);
    impl::apply(grad, phi);
}

/** Compute the divergence of a vector field
 *  \ingroup fvm
 */
template <typename FType>
std::unique_ptr<ScratchField> divergence(const FType& phi)
{
    AMREX_ALWAYS_ASSERT(phi.num_comp() >= AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(phi.num_comp() % AMREX_SPACEDIM == 0);
    const std::string gname = phi.name() + "_divergence";
    auto divphi =
        phi.repo().create_scratch_field(gname, phi.num_comp() / AMREX_SPACEDIM);
    divergence(*divphi, phi);
    return divphi;
}

} // namespace kynema_sgf::fvm

#endif /* DIVERGENCE_H */
