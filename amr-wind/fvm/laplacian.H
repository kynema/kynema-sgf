#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "amr-wind/fvm/fvm_utils.H"

namespace amr_wind::fvm {

/** Laplacian operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Laplacian
{
    Laplacian(FTypeOut& lphi, const FTypeIn& phi) : m_lapphi(lphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(m_lapphi.num_comp() == 1);
        AMREX_ALWAYS_ASSERT(m_phi.num_comp() == AMREX_SPACEDIM);
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& lapphi = m_lapphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                amrex::Real sp1 = Stencil::s00;
                amrex::Real s = Stencil::s01;
                amrex::Real sm1 = Stencil::s02;
                amrex::Real d2phidx2 =
                    (sp1 * phi(i + 1, j, k, 0) + s * phi(i, j, k, 0) +
                     sm1 * phi(i - 1, j, k, 0)) *
                    idx[0] * idx[0];
                sp1 = Stencil::s10;
                s = Stencil::s11;
                sm1 = Stencil::s12;
                amrex::Real d2phidy2 =
                    (sp1 * phi(i, j + 1, k, 1) + s * phi(i, j, k, 1) +
                     sm1 * phi(i, j - 1, k, 1)) *
                    idx[1] * idx[1];
                sp1 = Stencil::s20;
                s = Stencil::s21;
                sm1 = Stencil::s22;
                amrex::Real d2phidz2 =
                    (sp1 * phi(i, j, k + 1, 2) + s * phi(i, j, k, 2) +
                     sm1 * phi(i, j, k - 1, 2)) *
                    idx[2] * idx[2];
                lapphi(i, j, k) = d2phidx2 + d2phidy2 + d2phidz2;
            }
        });
    }

    FTypeOut& m_lapphi;
    const FTypeIn& m_phi;
};

/** Compute the laplacian of a given field
 *  \ingroup fvm
 *
 *  \param lapphi [inout] Field where the laplacian term is populated
 *  \param phi [in] Field whose gradient is computed
 */
template <typename FTypeIn, typename FTypeOut>
void laplacian(FTypeOut& lapphi, const FTypeIn& phi)
{
    BL_PROFILE("amr-wind::fvm::laplacian");
    Laplacian<FTypeIn, FTypeOut> lap(lapphi, phi);
    impl::apply(lap, phi);
}

/** Compute the laplacian of a given field and return as ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] Field whose gradient is computed
 */
template <typename FType>
std::unique_ptr<ScratchField> laplacian(const FType& phi)
{
    const std::string gname = phi.name() + "_laplacian";
    auto lapphi = phi.repo().create_scratch_field(gname, 1);
    laplacian(*lapphi, phi);
    return lapphi;
}

} // namespace amr_wind::fvm

#endif /* LAPLACIAN_H */
