#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "src/fvm/fvm_utils.H"

namespace kynema_sgf::fvm {

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
    void apply(const int lev) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();

        auto& out_mf = m_lapphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi = in_arrs[nbx];
                auto const& lapphi = out_arrs[nbx];
                const amrex::Real d2phidx2 =
                    (Stencil::s00 * phi(i + 1, j, k, 0) +
                     Stencil::s01 * phi(i, j, k, 0) +
                     Stencil::s02 * phi(i - 1, j, k, 0)) *
                    idx[0] * idx[0];
                const amrex::Real d2phidy2 =
                    (Stencil::s10 * phi(i, j + 1, k, 1) +
                     Stencil::s11 * phi(i, j, k, 1) +
                     Stencil::s12 * phi(i, j - 1, k, 1)) *
                    idx[1] * idx[1];
                const amrex::Real d2phidz2 =
                    (Stencil::s20 * phi(i, j, k + 1, 2) +
                     Stencil::s21 * phi(i, j, k, 2) +
                     Stencil::s22 * phi(i, j, k - 1, 2)) *
                    idx[2] * idx[2];
                lapphi(i, j, k) = d2phidx2 + d2phidy2 + d2phidz2;
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
    BL_PROFILE("kynema-sgf::fvm::laplacian");
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

} // namespace kynema_sgf::fvm

#endif /* LAPLACIAN_H */
