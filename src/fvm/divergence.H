#ifndef DIVERGENCE_H
#define DIVERGENCE_H

#include "src/fvm/fvm_utils.H"

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
    void apply(const int lev) const
    {
        const int ncomp = m_divphi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();

        auto& out_mf = m_divphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi_arr = in_arrs[nbx];
                auto const& divphi_arr = out_arrs[nbx];
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    divphi_arr(i, j, k, icomp) =
                        (Stencil::c00 *
                             phi_arr(
                                 i + 1, j, k, (icomp * AMREX_SPACEDIM) + 0) +
                         Stencil::c01 *
                             phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 0) +
                         Stencil::c02 *
                             phi_arr(
                                 i - 1, j, k, (icomp * AMREX_SPACEDIM) + 0)) *
                        idx[0];

                    divphi_arr(i, j, k, icomp) +=
                        (Stencil::c10 *
                             phi_arr(
                                 i, j + 1, k, (icomp * AMREX_SPACEDIM) + 1) +
                         Stencil::c11 *
                             phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 1) +
                         Stencil::c12 *
                             phi_arr(
                                 i, j - 1, k, (icomp * AMREX_SPACEDIM) + 1)) *
                        idx[1];

                    divphi_arr(i, j, k, icomp) +=
                        (Stencil::c20 *
                             phi_arr(
                                 i, j, k + 1, (icomp * AMREX_SPACEDIM) + 2) +
                         Stencil::c21 *
                             phi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 2) +
                         Stencil::c22 *
                             phi_arr(
                                 i, j, k - 1, (icomp * AMREX_SPACEDIM) + 2)) *
                        idx[2];
                }
            });
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
