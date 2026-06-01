#ifndef VORTICITYMAG_H
#define VORTICITYMAG_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::fvm {

/** VorticityMag operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct VorticityMag
{
    VorticityMag(FTypeOut& vortmagphi, const FTypeIn& phi)
        : m_vortmagphi(vortmagphi), m_phi(phi)
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

        auto& out_mf = m_vortmagphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi = in_arrs[nbx];
                auto const& vortmagphi = out_arrs[nbx];
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
                vortmagphi(i, j, k) = std::sqrt(
                    utils::powi(uy - vx, 2) + utils::powi(vz - wy, 2) +
                    utils::powi(wx - uz, 2));
            });
    }

    FTypeOut& m_vortmagphi;
    const FTypeIn& m_phi;
};

/** Compute the magnitude of vorticity
 *  \ingroup fvm
 *
 *  \param vortmagphi [out] Field where the vorticity magnitude is populated
 *  \param phi [in] The velocity vector field
 */
template <typename FTypeIn, typename FTypeOut>
void vorticity_mag(FTypeOut& vortmagphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::vorticity_mag");
    VorticityMag<FTypeIn, FTypeOut> vortmag(vortmagphi, phi);
    impl::apply(vortmag, phi);
}

/** Compute the magnitude of vorticity return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
std::unique_ptr<ScratchField> vorticity_mag(const FType& phi)
{
    const std::string gname = phi.name() + "_vorticity_mag";
    auto vortmagphi = phi.repo().create_scratch_field(gname, 1);
    vorticity_mag(*vortmagphi, phi);
    return vortmagphi;
}

} // namespace kynema_sgf::fvm

#endif /* VORTICITYMAG_H */
