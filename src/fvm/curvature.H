#ifndef CURVATURE_H
#define CURVATURE_H

#include "amr-wind/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::fvm {

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
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& curphi = m_curphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }

        amrex::ParallelFor(
            bx, m_phi.num_comp(),
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                amrex::Real cp1, c, cm1;
                amrex::Real sp1, s, sm1;

                sp1 = Stencil::s00;
                s = Stencil::s01;
                sm1 = Stencil::s02;
                const amrex::Real phixx =
                    (sp1 * phi(i + 1, j, k, n) + s * phi(i, j, k, n) +
                     sm1 * phi(i - 1, j, k, n)) *
                    idx[0] * idx[0];

                sp1 = Stencil::s10;
                s = Stencil::s11;
                sm1 = Stencil::s12;
                const amrex::Real phiyy =
                    (sp1 * phi(i, j + 1, k, n) + s * phi(i, j, k, n) +
                     sm1 * phi(i, j - 1, k, n)) *
                    idx[1] * idx[1];

                sp1 = Stencil::s20;
                s = Stencil::s21;
                sm1 = Stencil::s22;
                const amrex::Real phizz =
                    (sp1 * phi(i, j, k + 1, n) + s * phi(i, j, k, n) +
                     sm1 * phi(i, j, k - 1, n)) *
                    idx[2] * idx[2];

                cp1 = Stencil::c20;
                c = Stencil::c21;
                cm1 = Stencil::c22;
                const amrex::Real phiz =
                    (cp1 * phi(i, j, k + 1, n) + c * phi(i, j, k, n) +
                     cm1 * phi(i, j, k - 1, n)) *
                    idx[2];
                const amrex::Real phiz_ip1 =
                    (cp1 * phi(i + 1, j, k + 1, n) + c * phi(i + 1, j, k, n) +
                     cm1 * phi(i + 1, j, k - 1, n)) *
                    idx[2];
                const amrex::Real phiz_im1 =
                    (cp1 * phi(i - 1, j, k + 1, n) + c * phi(i - 1, j, k, n) +
                     cm1 * phi(i - 1, j, k - 1, n)) *
                    idx[2];
                const amrex::Real phiz_jp1 =
                    (cp1 * phi(i, j + 1, k + 1, n) + c * phi(i, j + 1, k, n) +
                     cm1 * phi(i, j + 1, k - 1, n)) *
                    idx[2];
                const amrex::Real phiz_jm1 =
                    (cp1 * phi(i, j - 1, k + 1, n) + c * phi(i, j - 1, k, n) +
                     cm1 * phi(i, j - 1, k - 1, n)) *
                    idx[2];

                cp1 = Stencil::c10;
                c = Stencil::c11;
                cm1 = Stencil::c12;
                const amrex::Real phiy =
                    (cp1 * phi(i, j + 1, k, n) + c * phi(i, j, k, n) +
                     cm1 * phi(i, j - 1, k, n)) *
                    idx[1];
                const amrex::Real phiy_ip1 =
                    (cp1 * phi(i + 1, j + 1, k, n) + c * phi(i + 1, j, k, n) +
                     cm1 * phi(i + 1, j - 1, k, n)) *
                    idx[1];
                const amrex::Real phiy_im1 =
                    (cp1 * phi(i - 1, j + 1, k, n) + c * phi(i - 1, j, k, n) +
                     cm1 * phi(i - 1, j - 1, k, n)) *
                    idx[1];
                const amrex::Real phiyz =
                    (cp1 * phiz_jp1 + c * phiz + cm1 * phiz_jm1) * idx[1];

                cp1 = Stencil::c00;
                c = Stencil::c01;
                cm1 = Stencil::c02;
                const amrex::Real phix =
                    (cp1 * phi(i + 1, j, k, n) + c * phi(i, j, k, n) +
                     cm1 * phi(i - 1, j, k, n)) *
                    idx[0];
                const amrex::Real phixy =
                    (cp1 * phiy_ip1 + c * phiy + cm1 * phiy_im1) * idx[0];
                const amrex::Real phixz =
                    (cp1 * phiz_ip1 + c * phiz + cm1 * phiz_im1) * idx[0];

                curphi(i, j, k, n) =
                    -((phix * phix * phiyy) - (2.0_rt * phix * phiy * phixy) +
                      (phiy * phiy * phixx) + (phix * phix * phizz) -
                      (2.0_rt * phix * phiz * phixz) + (phiz * phiz * phixx) +
                      (phiy * phiy * phizz) - (2.0_rt * phiy * phiz * phiyz) +
                      (phiz * phiz * phiyy)) /
                    std::pow(
                        (phix * phix) + (phiy * phiy) + (phiz * phiz), 1.5_rt);
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
    BL_PROFILE("amr-wind::fvm::curvature");
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

} // namespace amr_wind::fvm

#endif /* CURVATURE_H */
