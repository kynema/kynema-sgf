#ifndef CURVATURE_H
#define CURVATURE_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::fvm {

/** Curvature operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Curvature
{
    Curvature(FTypeOut& curphi, const FTypeIn& phi)
        : m_curphi(curphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(m_phi.num_comp() == m_curphi.num_comp());
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();
        const int ncomp = m_phi.num_comp();

        auto& out_mf = m_curphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi = in_arrs[nbx];
                auto const& curphi = out_arrs[nbx];
                for (int n = 0; n < ncomp; ++n) {
                    const amrex::Real phixx =
                        (Stencil::s00 * phi(i + 1, j, k, n) +
                         Stencil::s01 * phi(i, j, k, n) +
                         Stencil::s02 * phi(i - 1, j, k, n)) *
                        idx[0] * idx[0];
                    const amrex::Real phiyy =
                        (Stencil::s10 * phi(i, j + 1, k, n) +
                         Stencil::s11 * phi(i, j, k, n) +
                         Stencil::s12 * phi(i, j - 1, k, n)) *
                        idx[1] * idx[1];
                    const amrex::Real phizz =
                        (Stencil::s20 * phi(i, j, k + 1, n) +
                         Stencil::s21 * phi(i, j, k, n) +
                         Stencil::s22 * phi(i, j, k - 1, n)) *
                        idx[2] * idx[2];
                    const amrex::Real phiz =
                        (Stencil::c20 * phi(i, j, k + 1, n) +
                         Stencil::c21 * phi(i, j, k, n) +
                         Stencil::c22 * phi(i, j, k - 1, n)) *
                        idx[2];
                    const amrex::Real phiz_ip1 =
                        (Stencil::c20 * phi(i + 1, j, k + 1, n) +
                         Stencil::c21 * phi(i + 1, j, k, n) +
                         Stencil::c22 * phi(i + 1, j, k - 1, n)) *
                        idx[2];
                    const amrex::Real phiz_im1 =
                        (Stencil::c20 * phi(i - 1, j, k + 1, n) +
                         Stencil::c21 * phi(i - 1, j, k, n) +
                         Stencil::c22 * phi(i - 1, j, k - 1, n)) *
                        idx[2];
                    const amrex::Real phiz_jp1 =
                        (Stencil::c20 * phi(i, j + 1, k + 1, n) +
                         Stencil::c21 * phi(i, j + 1, k, n) +
                         Stencil::c22 * phi(i, j + 1, k - 1, n)) *
                        idx[2];
                    const amrex::Real phiz_jm1 =
                        (Stencil::c20 * phi(i, j - 1, k + 1, n) +
                         Stencil::c21 * phi(i, j - 1, k, n) +
                         Stencil::c22 * phi(i, j - 1, k - 1, n)) *
                        idx[2];
                    const amrex::Real phiy =
                        (Stencil::c10 * phi(i, j + 1, k, n) +
                         Stencil::c11 * phi(i, j, k, n) +
                         Stencil::c12 * phi(i, j - 1, k, n)) *
                        idx[1];
                    const amrex::Real phiy_ip1 =
                        (Stencil::c10 * phi(i + 1, j + 1, k, n) +
                         Stencil::c11 * phi(i + 1, j, k, n) +
                         Stencil::c12 * phi(i + 1, j - 1, k, n)) *
                        idx[1];
                    const amrex::Real phiy_im1 =
                        (Stencil::c10 * phi(i - 1, j + 1, k, n) +
                         Stencil::c11 * phi(i - 1, j, k, n) +
                         Stencil::c12 * phi(i - 1, j - 1, k, n)) *
                        idx[1];
                    const amrex::Real phiyz =
                        (Stencil::c10 * phiz_jp1 + Stencil::c11 * phiz +
                         Stencil::c12 * phiz_jm1) *
                        idx[1];
                    const amrex::Real phix =
                        (Stencil::c00 * phi(i + 1, j, k, n) +
                         Stencil::c01 * phi(i, j, k, n) +
                         Stencil::c02 * phi(i - 1, j, k, n)) *
                        idx[0];
                    const amrex::Real phixy =
                        (Stencil::c00 * phiy_ip1 + Stencil::c01 * phiy +
                         Stencil::c02 * phiy_im1) *
                        idx[0];
                    const amrex::Real phixz =
                        (Stencil::c00 * phiz_ip1 + Stencil::c01 * phiz +
                         Stencil::c02 * phiz_im1) *
                        idx[0];
                    curphi(i, j, k, n) =
                        -((phix * phix * phiyy) -
                          (2.0_rt * phix * phiy * phixy) +
                          (phiy * phiy * phixx) + (phix * phix * phizz) -
                          (2.0_rt * phix * phiz * phixz) +
                          (phiz * phiz * phixx) + (phiy * phiy * phizz) -
                          (2.0_rt * phiy * phiz * phiyz) +
                          (phiz * phiz * phiyy)) /
                        std::pow(
                            (phix * phix) + (phiy * phiy) + (phiz * phiz),
                            1.5_rt);
                }
            });
    }

    FTypeOut& m_curphi;
    const FTypeIn& m_phi;
};

/** Compute curvature
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
void curvature(FTypeOut& curphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::curvature");
    Curvature<FTypeIn, FTypeOut> cur(curphi, phi);
    impl::apply(cur, phi);
}

/** Compute curvature
 *  \ingroup fvm
 */
template <typename FType>
std::unique_ptr<ScratchField> curvature(const FType& phi)
{
    const std::string gname = phi.name() + "_curvature";
    auto curphi = phi.repo().create_scratch_field(gname, phi.num_comp());
    curvature(*curphi, phi);
    return curphi;
}

} // namespace kynema_sgf::fvm

#endif /* CURVATURE_H */
