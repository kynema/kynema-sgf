#ifndef FILTER_H
#define FILTER_H

#include "src/fvm/fvm_utils.H"
#include "AMReX_Array4.H"
#include "AMReX_Geometry.H"
#include "AMReX_MFIter.H"
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
    void apply(const int lev, const amrex::MFIter& mfi) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto& filterphi_arr = m_filterphi(lev).array(mfi);
        const auto& phi_arr = m_phi(lev).const_array(mfi);

        const auto& bx_in = mfi.tilebox();
        const auto& bx = Stencil::box(bx_in, geom);
        if (bx.isEmpty()) {
            return;
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            for (int icomp = 0; icomp < ncomp; icomp++) {
                amrex::Real fp1 = Stencil::f00;
                amrex::Real f = Stencil::f01;
                amrex::Real fm1 = Stencil::f02;
                const amrex::Real filx =
                    ((fp1 * phi_arr(i + 1, j, k, icomp)) +
                     (f * phi_arr(i, j, k, icomp)) +
                     (fm1 * phi_arr(i - 1, j, k, icomp)));

                fp1 = Stencil::f10;
                f = Stencil::f11;
                fm1 = Stencil::f12;
                const amrex::Real fily =
                    ((fp1 * phi_arr(i, j + 1, k, icomp)) +
                     (f * phi_arr(i, j, k, icomp)) +
                     (fm1 * phi_arr(i, j - 1, k, icomp)));

                fp1 = Stencil::f20;
                f = Stencil::f21;
                fm1 = Stencil::f22;
                const amrex::Real filz =
                    ((fp1 * phi_arr(i, j, k + 1, icomp)) +
                     (f * phi_arr(i, j, k, icomp)) +
                     (fm1 * phi_arr(i, j, k - 1, icomp)));
                filterphi_arr(i, j, k, icomp) =
                    1.0_rt / 3.0_rt * (filx + fily + filz);
            }
        });
    }

    void applyFused(const int lev) const
    {
        const int ncomp = m_phi.num_comp();
        const auto& geom = m_phi.repo().mesh().Geom(lev);
        const auto dlo = geom.Domain().smallEnd();
        const auto dhi = geom.Domain().bigEnd();
        const amrex::GpuArray<int, AMREX_SPACEDIM> periodic{
            geom.isPeriodic(0) ? 1 : 0,
            geom.isPeriodic(1) ? 1 : 0,
            geom.isPeriodic(2) ? 1 : 0};

        auto& out_mf = m_filterphi(lev);
        auto const& in_arrs = m_phi(lev).const_arrays();
        auto const& out_arrs = out_mf.arrays();

        amrex::ParallelFor(
            out_mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const auto sc = stencil::pick_stencil(
                    i, j, k, dlo, dhi, periodic);
                auto const& phi_arr = in_arrs[nbx];
                auto const& filterphi_arr = out_arrs[nbx];
                for (int icomp = 0; icomp < ncomp; icomp++) {
                    const amrex::Real filx =
                        sc.f00 * phi_arr(i + 1, j, k, icomp) +
                        sc.f01 * phi_arr(i, j, k, icomp) +
                        sc.f02 * phi_arr(i - 1, j, k, icomp);
                    const amrex::Real fily =
                        sc.f10 * phi_arr(i, j + 1, k, icomp) +
                        sc.f11 * phi_arr(i, j, k, icomp) +
                        sc.f12 * phi_arr(i, j - 1, k, icomp);
                    const amrex::Real filz =
                        sc.f20 * phi_arr(i, j, k + 1, icomp) +
                        sc.f21 * phi_arr(i, j, k, icomp) +
                        sc.f22 * phi_arr(i, j, k - 1, icomp);
                    filterphi_arr(i, j, k, icomp) =
                        1.0_rt / 3.0_rt * (filx + fily + filz);
                }
            });
        amrex::Gpu::streamSynchronize();
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
