#ifndef QCRITERION_H
#define QCRITERION_H

#include "amr-wind/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::fvm {

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

            const amrex::Real S2 = std::pow(ux, 2.0_rt) + std::pow(vy, 2.0_rt) +
                                   std::pow(wz, 2.0_rt) +
                                   (0.5_rt * std::pow(uy + vx, 2.0_rt)) +
                                   (0.5_rt * std::pow(vz + wy, 2.0_rt)) +
                                   (0.5_rt * std::pow(wx + uz, 2.0_rt));

            const amrex::Real W2 = (0.5_rt * std::pow(uy - vx, 2.0_rt)) +
                                   (0.5_rt * std::pow(vz - wy, 2.0_rt)) +
                                   (0.5_rt * std::pow(wx - uz, 2.0_rt));
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
    BL_PROFILE("amr-wind::fvm::q_criterion");
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

} // namespace amr_wind::fvm

#endif /* QCRITERION_H */
