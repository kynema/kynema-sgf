#ifndef FILTER_H
#define FILTER_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::fvm {

/** filter operator
 *  \ingroup fvm
 */
template <typename FTypeIn, typename FTypeOut>
struct Filter
{
    /**
     *  \param filterphi The filter field \f$\widetilde{ \mathbf{\phi}}\f$
     *  \param phi The input file \f$\phi\f$
     */
    Filter(FTypeOut& filterphi, const FTypeIn& phi)
        : m_filterphi(filterphi), m_phi(phi)
    {
        AMREX_ALWAYS_ASSERT(filterphi.num_comp() == phi.num_comp());
        AMREX_ALWAYS_ASSERT(m_phi.num_grow() > amrex::IntVect(0));
    }

    template <typename Stencil>
    void apply(const int lev) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();

        auto& out_mf = m_filterphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if (!Stencil::applies_to(i, j, k, dlo, dhi)) {
                    return;
                }
                auto const& phi_arr = in_arrs[nbx];
                auto const& filterphi_arr = out_arrs[nbx];
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    const amrex::Real filx =
                        Stencil::f00 * phi_arr(i + 1, j, k, icomp) +
                        Stencil::f01 * phi_arr(i, j, k, icomp) +
                        Stencil::f02 * phi_arr(i - 1, j, k, icomp);
                    const amrex::Real fily =
                        Stencil::f10 * phi_arr(i, j + 1, k, icomp) +
                        Stencil::f11 * phi_arr(i, j, k, icomp) +
                        Stencil::f12 * phi_arr(i, j - 1, k, icomp);
                    const amrex::Real filz =
                        Stencil::f20 * phi_arr(i, j, k + 1, icomp) +
                        Stencil::f21 * phi_arr(i, j, k, icomp) +
                        Stencil::f22 * phi_arr(i, j, k - 1, icomp);
                    filterphi_arr(i, j, k, icomp) =
                        1.0_rt / 3.0_rt * (filx + fily + filz);
                }
            });
    }

    FTypeOut& m_filterphi;
    const FTypeIn& m_phi;
};

/** Compute the filtered value of a given field
 *  \ingroup fvm
 *
 *  \param filterphi [inout] Field where the filtered term is populated
 *  \param phi [in] Field whose filter is computed
 */
template <typename FTypeIn, typename FTypeOut>
void filter(FTypeOut& filterphi, const FTypeIn& phi)
{
    BL_PROFILE("kynema-sgf::fvm::filter");
    Filter<FTypeIn, FTypeOut> filtering(filterphi, phi);
    impl::apply(filtering, phi);
}

/** Compute the filter of a given field and return it as ScratchField
 *  \ingroup fvm
 *
 *  \param phi [in] Field whose filter is computed
 */
template <typename FType>
std::unique_ptr<ScratchField> filter(const FType& phi)
{
    const std::string gname = phi.name() + "_filter";
    auto filterphi = phi.repo().create_scratch_field(gname, phi.num_comp());
    filter(*filterphi, phi);
    return filterphi;
}

} // namespace kynema_sgf::fvm

#endif /* FILTER_H */
