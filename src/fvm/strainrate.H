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
                strphi(i, j, k) = std::sqrt(
                    (2.0_rt * utils::powi(ux, 2)) +
                    (2.0_rt * utils::powi(vy, 2)) +
                    (2.0_rt * utils::powi(wz, 2)) + utils::powi(uy + vx, 2) +
                    utils::powi(vz + wy, 2) + utils::powi(wx + uz, 2));
            });
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
