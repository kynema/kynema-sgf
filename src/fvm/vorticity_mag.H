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
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& vortmagphi = m_vortmagphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real cp1, c, cm1, uy, uz, vx, vz, wx, wy;
            cp1 = Stencil::c00;
            c = Stencil::c01;
            cm1 = Stencil::c02;

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

            vortmagphi(i, j, k) = std::sqrt(
                utils::powi(uy - vx, 2) + utils::powi(vz - wy, 2) +
                utils::powi(wx - uz, 2));
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

        auto& out_mf = m_vortmagphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const auto sc = stencil::pick_stencil(
                    i, j, k, dlo, dhi, periodic);
                auto const& phi = in_arrs[nbx];
                auto const& vortmagphi = out_arrs[nbx];
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
                vortmagphi(i, j, k) = std::sqrt(
                    utils::powi(uy - vx, 2) + utils::powi(vz - wy, 2) +
                    utils::powi(wx - uz, 2));
            });
        amrex::Gpu::streamSynchronize();
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
