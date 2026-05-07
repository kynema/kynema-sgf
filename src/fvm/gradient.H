#ifndef GRADIENT_H
#define GRADIENT_H

#include "src/fvm/fvm_utils.H"

namespace kynema_sgf::fvm {

/** Gradient operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Gradient
{
    /**
     *  \param gradphi The gradient field \f$\nabla \mathbf{\phi}\f$
     *  \param phi The input file \f$\phi\f$
     */
    Gradient(FTypeOut& gradphi, const FTypeIn& phi)
        : m_gradphi(gradphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(
            gradphi.num_comp() == AMREX_SPACEDIM * phi.num_comp());
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();

        auto& out_mf = m_gradphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi_arr = in_arrs[nbx];
                auto const& gradphi_arr = out_arrs[nbx];
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    gradphi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 0) =
                        (Stencil::c00 * phi_arr(i + 1, j, k, icomp) +
                         Stencil::c01 * phi_arr(i, j, k, icomp) +
                         Stencil::c02 * phi_arr(i - 1, j, k, icomp)) *
                        idx[0];

                    gradphi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 1) =
                        (Stencil::c10 * phi_arr(i, j + 1, k, icomp) +
                         Stencil::c11 * phi_arr(i, j, k, icomp) +
                         Stencil::c12 * phi_arr(i, j - 1, k, icomp)) *
                        idx[1];

                    gradphi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 2) =
                        (Stencil::c20 * phi_arr(i, j, k + 1, icomp) +
                         Stencil::c21 * phi_arr(i, j, k, icomp) +
                         Stencil::c22 * phi_arr(i, j, k - 1, icomp)) *
                        idx[2];
                }
            });
    }

    FTypeOut& m_gradphi;
    const FTypeIn& m_phi;
};

/** Compute the gradient of a given field
 *  \ingroup fvm
 *
 *  \param gradphi [inout] Field where the gradient term is populated
 *  \param phi [in] Field whose gradient is computed
 */
template <typename FTypeIn, typename FTypeOut>
void gradient(FTypeOut& gradphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::gradient");
    Gradient<FTypeIn, FTypeOut> grad(gradphi, phi);
    impl::apply(grad, phi);
}

/** Compute the gradient of a given field and return it as ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] Field whose gradient is computed
 */
template <typename FType>
std::unique_ptr<ScratchField> gradient(const FType& phi)
{
    const std::string gname = phi.name() + "_gradient";
    auto gradphi =
        phi.repo().create_scratch_field(gname, phi.num_comp() * AMREX_SPACEDIM);
    gradient(*gradphi, phi);
    return gradphi;
}

} // namespace kynema_sgf::fvm

#endif /* GRADIENT_H */
