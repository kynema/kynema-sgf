#ifndef ADVOP_MOL_H
#define ADVOP_MOL_H

#include <type_traits>

#include "amr-wind/convection/MOL.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"

#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::pde {

/** MOL advection operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct AdvectionOp<
    PDE,
    fvm::MOL,
    typename std::enable_if_t<std::is_base_of_v<ScalarTransport, PDE>>>
{
    // cppcheck-suppress uninitMemberVar
    AdvectionOp(
        CFDSim& /* sim */,
        PDEFields& fields_in,
        bool /* has_overset */,
        bool /* variable density */,
        bool /* mesh mapping */,
        bool /* is_anelastic */)
        : fields(fields_in)
        , density(fields_in.repo.get_field("density"))
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
    {}

    void preadvect(
        const FieldState /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/)
    {}

    void operator()(const FieldState fstate, const amrex::Real /*unused*/)
    {
        static_assert(
            PDE::ndim == 1, "Invalid number of components for scalar");

        const auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        // cppcheck-suppress constVariableReference
        auto& conv_term = fields.conv_term.state(fstate);
        const auto& dof_field = fields.field.state(fstate);
        const auto& den = density.state(fstate);

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MFItInfo mfi_info;
            // if (amrex::Gpu::notInLaunchRegion())
            // mfi_info.EnableTiling(amrex::IntVect(1024,16,16)).SetDynamic(true);
            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(density(lev), mfi_info); mfi.isValid();
                 ++mfi) {

                amrex::Box const& bx = mfi.tilebox();
                auto rho_arr = den(lev).const_array(mfi);
                auto tra_arr = dof_field(lev).const_array(mfi);
                amrex::FArrayBox rhotracfab;
                amrex::Array4<amrex::Real> rhotrac;

                if (PDE::multiply_rho) {

                    amrex::Box rhotrac_box =
                        amrex::grow(bx, fvm::MOL::nghost_state);
                    rhotracfab.resize(
                        rhotrac_box, PDE::ndim, amrex::The_Async_Arena());
                    rhotrac = rhotracfab.array();

                    amrex::ParallelFor(
                        rhotrac_box, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                            rhotrac(i, j, k, n) =
                                rho_arr(i, j, k) * tra_arr(i, j, k, n);
                        });
                }

                {
                    const int nmaxcomp = PDE::ndim;

                    amrex::Box tmpbox = amrex::surroundingNodes(bx);
                    const int tmpcomp = nmaxcomp * AMREX_SPACEDIM;

                    amrex::FArrayBox tmpfab(
                        tmpbox, tmpcomp, amrex::The_Async_Arena());

                    amrex::Array4<amrex::Real> fx = tmpfab.array(0);
                    amrex::Array4<amrex::Real> fy = tmpfab.array(nmaxcomp);
                    amrex::Array4<amrex::Real> fz = tmpfab.array(nmaxcomp * 2);

                    mol::compute_convective_fluxes(
                        lev, bx, PDE::ndim, fx, fy, fz,
                        (PDE::multiply_rho ? rhotrac : tra_arr),
                        u_mac(lev).const_array(mfi),
                        v_mac(lev).const_array(mfi),
                        w_mac(lev).const_array(mfi), dof_field.bcrec().data(),
                        dof_field.bcrec_device().data(), geom);

                    mol::compute_convective_rate(
                        bx, PDE::ndim, conv_term(lev).array(mfi), fx, fy, fz,
                        geom[lev].InvCellSizeArray());
                }
            }
        }
    }

    PDEFields& fields;
    Field& density;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
};

} // namespace amr_wind::pde

#endif /* ADVOP_MOL_H */
