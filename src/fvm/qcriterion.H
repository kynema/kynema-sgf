#ifndef QCRITERION_H
#define QCRITERION_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::fvm {

/** Qcriterion operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Qcriterion
{
    Qcriterion(
        FTypeOut& qcritphi, const FTypeIn& phi, const bool nondim = false)
        : m_qcritphi(qcritphi), m_phi(phi), m_nondim(nondim)
    {
        AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == m_phi.num_comp());
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto& qcritphi = m_qcritphi(lev).array(mfi);
        const auto& phi = m_phi(lev).const_array(mfi);
        const bool nondim = m_nondim;

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

            const amrex::Real S2 = utils::powi(ux, 2) + utils::powi(vy, 2) +
                                   utils::powi(wz, 2) +
                                   (0.5_rt * utils::powi(uy + vx, 2)) +
                                   (0.5_rt * utils::powi(vz + wy, 2)) +
                                   (0.5_rt * utils::powi(wx + uz, 2));

            const amrex::Real W2 = (0.5_rt * utils::powi(uy - vx, 2)) +
                                   (0.5_rt * utils::powi(vz - wy, 2)) +
                                   (0.5_rt * utils::powi(wx - uz, 2));
            if (nondim) {
                qcritphi(i, j, k) =
                    0.5_rt *
                    (W2 / amrex::max<amrex::Real>(
                              std::numeric_limits<amrex::Real>::epsilon() *
                                  1.0e2_rt,
                              S2) -
                     1.0_rt);
            } else {
                qcritphi(i, j, k) = 0.5_rt * (W2 - S2);
            }
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
        const bool nondim = m_nondim;

        auto& out_mf = m_qcritphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const auto sc = stencil::pick_stencil(
                    i, j, k, dlo, dhi, periodic);
                auto const& phi = in_arrs[nbx];
                auto const& qcritphi = out_arrs[nbx];
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
                const amrex::Real S2 =
                    utils::powi(ux, 2) + utils::powi(vy, 2) +
                    utils::powi(wz, 2) + (0.5_rt * utils::powi(uy + vx, 2)) +
                    (0.5_rt * utils::powi(vz + wy, 2)) +
                    (0.5_rt * utils::powi(wx + uz, 2));
                const amrex::Real W2 =
                    (0.5_rt * utils::powi(uy - vx, 2)) +
                    (0.5_rt * utils::powi(vz - wy, 2)) +
                    (0.5_rt * utils::powi(wx - uz, 2));
                if (nondim) {
                    qcritphi(i, j, k) =
                        0.5_rt *
                        (W2 / amrex::max<amrex::Real>(
                                  std::numeric_limits<amrex::Real>::epsilon() *
                                      1.0e2_rt,
                                  S2) -
                         1.0_rt);
                } else {
                    qcritphi(i, j, k) = 0.5_rt * (W2 - S2);
                }
            });
        amrex::Gpu::streamSynchronize();
    }

    FTypeOut& m_qcritphi;
    const FTypeIn& m_phi;
    bool m_nondim;
};

/** Compute the q-criterion
 *  \ingroup fvm
 *
 *  \param qcritphi [out] Field where the Q-criterion is populated
 *  \param phi [in] The velocity vector field
 *  \param nondim [in] Flag for non-dimensional value
 */
template <typename FTypeIn, typename FTypeOut>
void q_criterion(
    FTypeOut& qcritphi, const FTypeIn& phi, const bool nondim = false)
{
    BL_PROFILE("kynema-sgf::fvm::q_criterion");
    Qcriterion<FTypeIn, FTypeOut> qcriterion(qcritphi, phi, nondim);
    impl::apply(qcriterion, phi);
}

/** Compute the Q-criterion return as a ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] The velocity vector field
 */
template <typename FType>
std::unique_ptr<ScratchField> q_criterion(const FType& phi)
{
    const std::string gname = phi.name() + "_q_criterion";
    auto qcritphi = phi.repo().create_scratch_field(gname, 1);
    q_criterion(*qcritphi, phi);
    return qcritphi;
}

} // namespace kynema_sgf::fvm

#endif /* QCRITERION_H */
