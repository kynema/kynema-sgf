#ifndef GRADIENT_H
#define GRADIENT_H

#include "amr-wind/fvm/fvm_utils.H"
#include "AMReX_Array4.H"
#include "AMReX_Geometry.H"
#include "AMReX_MFIter.H"

namespace amr_wind::fvm {

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
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& gradphi_arr = m_gradphi(lev).array(mfi);
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
                gradphi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 0) =
                    (cp1 * phi_arr(i + 1, j, k, icomp) +
                     c * phi_arr(i, j, k, icomp) +
                     cm1 * phi_arr(i - 1, j, k, icomp)) *
                    idx[0];

                cp1 = Stencil::c10;
                c = Stencil::c11;
                cm1 = Stencil::c12;
                gradphi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 1) =
                    (cp1 * phi_arr(i, j + 1, k, icomp) +
                     c * phi_arr(i, j, k, icomp) +
                     cm1 * phi_arr(i, j - 1, k, icomp)) *
                    idx[1];

                cp1 = Stencil::c20;
                c = Stencil::c21;
                cm1 = Stencil::c22;
                gradphi_arr(i, j, k, (icomp * AMREX_SPACEDIM) + 2) =
                    (cp1 * phi_arr(i, j, k + 1, icomp) +
                     c * phi_arr(i, j, k, icomp) +
                     cm1 * phi_arr(i, j, k - 1, icomp)) *
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
    BL_PROFILE("amr-wind::fvm::gradient");
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

} // namespace amr_wind::fvm

#endif /* GRADIENT_H */
