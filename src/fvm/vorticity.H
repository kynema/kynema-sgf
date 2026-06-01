#ifndef VORTICITY_H
#define VORTICITY_H

#include "src/fvm/fvm_utils.H"

namespace kynema_sgf::fvm {

/** Vorticity operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Vorticity
{
    Vorticity(FTypeOut& vortphi, const FTypeIn& phi)
        : m_vort(vortphi), m_phi(phi)
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

        auto& out_mf = m_vort(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi = in_arrs[nbx];
                auto const& vort = out_arrs[nbx];
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
                vort(i, j, k, 0) = wy - vz;
                vort(i, j, k, 1) = uz - wx;
                vort(i, j, k, 2) = vx - uy;
            });
    }

    FTypeOut& m_vort;
    const FTypeIn& m_phi;
};

/** Compute the magnitude of vorticity
 *  \ingroup fvm
 *
 *  \param vortphi [out] Field where the vorticity magnitude is populated
 *  \param phi [in] The velocity vector field
 */
template <typename FTypeIn, typename FTypeOut>
void vorticity(FTypeOut& vortphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::vorticity");
    Vorticity<FTypeIn, FTypeOut> vort(vortphi, phi);
    impl::apply(vort, phi);
}

/** Compute the magnitude of vorticity return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
std::unique_ptr<ScratchField> vorticity(const FType& phi)
{
    const std::string gname = phi.name() + "_vorticity";
    auto vortphi = phi.repo().create_scratch_field(gname, AMREX_SPACEDIM);
    vorticity(*vortphi, phi);
    return vortphi;
}

} // namespace kynema_sgf::fvm

#endif /* VORTICITY_H */
