#ifndef NONLINEARSUM_H
#define NONLINEARSUM_H

#include "src/fvm/fvm_utils.H"

namespace kynema_sgf::fvm {

/** Non-linear Term for Kosovic Model
 *  \ingroup fvm
 *  \defgroup NonLinearSGS Non-linear Term for Kosovic Model
 */
template <typename FTypeIn, typename FTypeOut>
struct NonLinearSum
{
    NonLinearSum(FTypeOut& strphi, const FTypeIn& phi)
        : m_strphi(strphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == m_phi.num_comp());
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();

        auto& out_mf = m_strphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi = in_arrs[nbx];
                auto const& strphi = out_arrs[nbx];
                const amrex::Real ux = (Stencil::c00 * phi(i + 1, j, k, 0) +
                                        Stencil::c01 * phi(i, j, k, 0) +
                                        Stencil::c02 * phi(i - 1, j, k, 0)) *
                                       idx[0];
                const amrex::Real vx = (Stencil::c00 * phi(i + 1, j, k, 1) +
                                        Stencil::c01 * phi(i, j, k, 1) +
                                        Stencil::c02 * phi(i - 1, j, k, 1)) *
                                       idx[0];
                const amrex::Real wx = (Stencil::c00 * phi(i + 1, j, k, 2) +
                                        Stencil::c01 * phi(i, j, k, 2) +
                                        Stencil::c02 * phi(i - 1, j, k, 2)) *
                                       idx[0];
                const amrex::Real uy = (Stencil::c10 * phi(i, j + 1, k, 0) +
                                        Stencil::c11 * phi(i, j, k, 0) +
                                        Stencil::c12 * phi(i, j - 1, k, 0)) *
                                       idx[1];
                const amrex::Real vy = (Stencil::c10 * phi(i, j + 1, k, 1) +
                                        Stencil::c11 * phi(i, j, k, 1) +
                                        Stencil::c12 * phi(i, j - 1, k, 1)) *
                                       idx[1];
                const amrex::Real wy = (Stencil::c10 * phi(i, j + 1, k, 2) +
                                        Stencil::c11 * phi(i, j, k, 2) +
                                        Stencil::c12 * phi(i, j - 1, k, 2)) *
                                       idx[1];
                const amrex::Real uz = (Stencil::c20 * phi(i, j, k + 1, 0) +
                                        Stencil::c21 * phi(i, j, k, 0) +
                                        Stencil::c22 * phi(i, j, k - 1, 0)) *
                                       idx[2];
                const amrex::Real vz = (Stencil::c20 * phi(i, j, k + 1, 1) +
                                        Stencil::c21 * phi(i, j, k, 1) +
                                        Stencil::c22 * phi(i, j, k - 1, 1)) *
                                       idx[2];
                const amrex::Real wz = (Stencil::c20 * phi(i, j, k + 1, 2) +
                                        Stencil::c21 * phi(i, j, k, 2) +
                                        Stencil::c22 * phi(i, j, k - 1, 2)) *
                                       idx[2];
                // N11
                strphi(i, j, k, 0) = (ux * ux + uy * vx + uz * wx) -
                                     (ux * ux + uy * uy + uz * uz) +
                                     (3 * (ux * ux + vx * vx + wx * wx)) +
                                     (ux * ux + vx * uy + wx * uz);
                // N12
                strphi(i, j, k, 1) = (ux * uy + uy * vy + uz * wy) -
                                     (ux * vx + uy * vy + uz * vz) +
                                     (3 * (ux * uy + vx * vy + wx * wy)) +
                                     (ux * vx + vx * vy + wx * vz);
                // N13
                strphi(i, j, k, 2) = (ux * uz + uy * vz + uz * wz) -
                                     (ux * wx + uy * wy + uz * wz) +
                                     (3 * (ux * uz + vx * vz + wx * wz)) +
                                     (ux * wx + vx * wy + wx * wz);
                // N21
                strphi(i, j, k, 3) = (vx * ux + vy * vx + vz * wx) -
                                     (vx * ux + vy * uy + vz * uz) +
                                     (3 * (uy * ux + vy * vx + wy * wx)) +
                                     (uy * ux + vy * uy + wy * uz);
                // N22
                strphi(i, j, k, 4) = (vx * uy + vy * vy + vz * wy) -
                                     (vx * vx + vy * vy + vz * vz) +
                                     (3 * (uy * vx + vy * vy + wy * wy)) +
                                     (uy * vx + vy * vy + wy * vz);
                // N23
                strphi(i, j, k, 5) = (vx * uz + vy * vz + vz * wz) -
                                     (vx * wx + vy * wy + vz * wz) +
                                     (3 * (uy * wx + vy * vz + wy * wz)) +
                                     (uy * wx + vy * wy + wy * wz);
                // N31
                strphi(i, j, k, 6) = (wx * ux + wy * vx + wz * wx) -
                                     (wx * ux + wy * uy + wz * uz) +
                                     (3 * (uz * ux + vz * vx + wz * wx)) +
                                     (uz * ux + vz * uy + wz * uz);
                // N32
                strphi(i, j, k, 7) = (wx * uy + wy * vy + wz * wy) -
                                     (wx * vx + wy * vy + wz * vz) +
                                     (3 * (uz * uy + vz * vy + wz * wy)) +
                                     (uz * vx + vz * vy + wz * vz);
                // N33
                strphi(i, j, k, 8) = (wx * uz + wy * vz + wz * wz) -
                                     (wx * wx + wy * wy + wz * wz) +
                                     (3 * (uz * uz + vz * vz + wz * wz)) +
                                     (uz * wx + vz * wy + wz * wz);
            });
    }

    FTypeOut& m_strphi;
    const FTypeIn& m_phi;
};

/** Compute the non-linear term for the Kosovic Sub-grid scale model
 *  \ingroup fvm
 *
 *  \param strphi [out] Field where the non-linear terms is populated
 *  \param phi [in] The stress tensor  field
 */
template <typename FTypeIn, typename FTypeOut>
void nonlinearsum(FTypeOut& strphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::nonlinearsum");
    NonLinearSum<FTypeIn, FTypeOut> str(strphi, phi);
    impl::apply(str, phi);
}

/** Compute the non-linear term as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The stress tensor  field
 */
template <typename FType>
std::unique_ptr<ScratchField> nonlinearsum(const FType& phi)
{
    const std::string gname = phi.name() + "_nonlinearsum";
    auto strphi = phi.repo().create_scratch_field(gname, 9);
    nonlinearsum(*strphi, phi);
    return strphi;
}

} // namespace kynema_sgf::fvm

#endif /* NONLINEARSUM_H */
