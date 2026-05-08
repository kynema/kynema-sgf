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
    void apply(const int lev) const
    {
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& idx = geom.InvCellSizeArray();
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();
        const bool nondim = m_nondim;

        auto& out_mf = m_qcritphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi = in_arrs[nbx];
                auto const& qcritphi = out_arrs[nbx];
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
