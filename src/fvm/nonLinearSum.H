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
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& strphi = m_strphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real cp1, c, cm1, ux, uy, uz, vx, vy, vz, wx, wy, wz;
            cp1 = Stencil::c00;
            c = Stencil::c01;
            cm1 = Stencil::c02;

            ux = (cp1 * phi(i + 1, j, k, 0) + c * phi(i, j, k, 0) +
                  cm1 * phi(i - 1, j, k, 0)) *
                 idx[0];
            vx = (cp1 * phi(i + 1, j, k, 1) + c * phi(i, j, k, 1) +
                  cm1 * phi(i - 1, j, k, 1)) *
                 idx[0];
            wx = (cp1 * phi(i + 1, j, k, 2) + c * phi(i, j, k, 2) +
                  cm1 * phi(i - 1, j, k, 2)) *
                 idx[0];

            cp1 = Stencil::c10;
            c = Stencil::c11;
            cm1 = Stencil::c12;

            uy = (cp1 * phi(i, j + 1, k, 0) + c * phi(i, j, k, 0) +
                  cm1 * phi(i, j - 1, k, 0)) *
                 idx[1];
            vy = (cp1 * phi(i, j + 1, k, 1) + c * phi(i, j, k, 1) +
                  cm1 * phi(i, j - 1, k, 1)) *
                 idx[1];
            wy = (cp1 * phi(i, j + 1, k, 2) + c * phi(i, j, k, 2) +
                  cm1 * phi(i, j - 1, k, 2)) *
                 idx[1];

            cp1 = Stencil::c20;
            c = Stencil::c21;
            cm1 = Stencil::c22;

            uz = (cp1 * phi(i, j, k + 1, 0) + c * phi(i, j, k, 0) +
                  cm1 * phi(i, j, k - 1, 0)) *
                 idx[2];
            vz = (cp1 * phi(i, j, k + 1, 1) + c * phi(i, j, k, 1) +
                  cm1 * phi(i, j, k - 1, 1)) *
                 idx[2];
            wz = (cp1 * phi(i, j, k + 1, 2) + c * phi(i, j, k, 2) +
                  cm1 * phi(i, j, k - 1, 2)) *
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

    void applyFused(const int lev) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();
        const amrex::GpuArray<int, AMREX_SPACEDIM> periodic{
            geom.isPeriodic(0) ? 1 : 0,
            geom.isPeriodic(1) ? 1 : 0,
            geom.isPeriodic(2) ? 1 : 0};

        auto& out_mf = m_strphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const auto sc = stencil::pick_stencil(
                    i, j, k, dlo, dhi, periodic);
                auto const& phi = in_arrs[nbx];
                auto const& strphi = out_arrs[nbx];
                const amrex::Real ux =
                    (sc.c00 * phi(i + 1, j, k, 0) + sc.c01 * phi(i, j, k, 0) +
                     sc.c02 * phi(i - 1, j, k, 0)) *
                    idx[0];
                const amrex::Real vx =
                    (sc.c00 * phi(i + 1, j, k, 1) + sc.c01 * phi(i, j, k, 1) +
                     sc.c02 * phi(i - 1, j, k, 1)) *
                    idx[0];
                const amrex::Real wx =
                    (sc.c00 * phi(i + 1, j, k, 2) + sc.c01 * phi(i, j, k, 2) +
                     sc.c02 * phi(i - 1, j, k, 2)) *
                    idx[0];
                const amrex::Real uy =
                    (sc.c10 * phi(i, j + 1, k, 0) + sc.c11 * phi(i, j, k, 0) +
                     sc.c12 * phi(i, j - 1, k, 0)) *
                    idx[1];
                const amrex::Real vy =
                    (sc.c10 * phi(i, j + 1, k, 1) + sc.c11 * phi(i, j, k, 1) +
                     sc.c12 * phi(i, j - 1, k, 1)) *
                    idx[1];
                const amrex::Real wy =
                    (sc.c10 * phi(i, j + 1, k, 2) + sc.c11 * phi(i, j, k, 2) +
                     sc.c12 * phi(i, j - 1, k, 2)) *
                    idx[1];
                const amrex::Real uz =
                    (sc.c20 * phi(i, j, k + 1, 0) + sc.c21 * phi(i, j, k, 0) +
                     sc.c22 * phi(i, j, k - 1, 0)) *
                    idx[2];
                const amrex::Real vz =
                    (sc.c20 * phi(i, j, k + 1, 1) + sc.c21 * phi(i, j, k, 1) +
                     sc.c22 * phi(i, j, k - 1, 1)) *
                    idx[2];
                const amrex::Real wz =
                    (sc.c20 * phi(i, j, k + 1, 2) + sc.c21 * phi(i, j, k, 2) +
                     sc.c22 * phi(i, j, k - 1, 2)) *
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
        amrex::Gpu::streamSynchronize();
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
