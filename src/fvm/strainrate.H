#ifndef STRAINRATE_H
#define STRAINRATE_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::fvm {

/** Strain rate operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct StrainRate
{
    StrainRate(FTypeOut& strphi, const FTypeIn& phi)
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

            strphi(i, j, k) = std::sqrt(
                (2.0_rt * utils::powi(ux, 2)) + (2.0_rt * utils::powi(vy, 2)) +
                (2.0_rt * utils::powi(wz, 2)) + utils::powi(uy + vx, 2) +
                utils::powi(vz + wy, 2) + utils::powi(wx + uz, 2));
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
                strphi(i, j, k) = std::sqrt(
                    (2.0_rt * utils::powi(ux, 2)) +
                    (2.0_rt * utils::powi(vy, 2)) +
                    (2.0_rt * utils::powi(wz, 2)) + utils::powi(uy + vx, 2) +
                    utils::powi(vz + wy, 2) + utils::powi(wx + uz, 2));
            });
        amrex::Gpu::streamSynchronize();
    }

    FTypeOut& m_strphi;
    const FTypeIn& m_phi;
};

/** Compute the magnitude of strain rate
 *  \ingroup fvm
 *
 *  \param strphi [out] Field where the strain rate magnitude is populated
 *  \param phi [in] The velocity vector field
 */
template <typename FTypeIn, typename FTypeOut>
void strainrate(FTypeOut& strphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::strainrate");
    StrainRate<FTypeIn, FTypeOut> str(strphi, phi);
    impl::apply(str, phi);
}

/** Compute the magnitude of strain rate return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
std::unique_ptr<ScratchField> strainrate(const FType& phi)
{
    const std::string gname = phi.name() + "_strainrate";
    auto strphi = phi.repo().create_scratch_field(gname, 1);
    strainrate(*strphi, phi);
    return strphi;
}

} // namespace kynema_sgf::fvm

#endif /* STRAINRATE_H */
