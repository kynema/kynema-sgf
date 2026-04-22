#ifndef DENSITY_OPS_H
#define DENSITY_OPS_H

#include "amr-wind/equation_systems/density/density.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

template <typename Scheme>
struct ComputeRHSOp<Density, Scheme>
{
    explicit ComputeRHSOp(PDEFields& fields_in) : fields(fields_in) {}

    void predictor_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real dt,
        bool /*mesh_mapping*/)
    {
        // Field states for diffusion and advection terms. In Godunov scheme
        // these terms only have one state.
        auto fstate = std::is_same_v<Scheme, fvm::Godunov> ? FieldState::New
                                                           : FieldState::Old;

        const int nlevels = fields.repo.num_active_levels();
        auto& field = fields.field;
        const auto& field_old = field.state(FieldState::Old);
        const auto& conv_term = fields.conv_term.state(fstate);

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& rho_arrs = field(lev).arrays();
            const auto& rho_o_arrs = field_old(lev).const_arrays();
            const auto& ddt_o_arrs = conv_term(lev).const_arrays();

            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    rho_arrs[nbx](i, j, k) = rho_o_arrs[nbx](i, j, k) +
                                             (dt * ddt_o_arrs[nbx](i, j, k));
                    // Defer n+1/2 update to the time-stepping algorithm
                    // rho_nph(i, j, k) = 0.5_rt * (rho(i, j, k) + rho_o(i, j,
                    // k));
                });
        }
        amrex::Gpu::streamSynchronize();
    }

    void corrector_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real dt,
        bool /*mesh_mapping*/)
    {
        const int nlevels = fields.repo.num_active_levels();
        auto& field = fields.field;
        const auto& field_old = field.state(FieldState::Old);
        const auto& conv_term = fields.conv_term;
        const auto& conv_term_old = fields.conv_term.state(FieldState::Old);

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& rho_arrs = field(lev).arrays();
            const auto& rho_o_arrs = field_old(lev).const_arrays();
            const auto& ddt_arrs = conv_term(lev).const_arrays();
            const auto& ddt_o_arrs = conv_term_old(lev).const_arrays();

            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    rho_arrs[nbx](i, j, k) =
                        rho_o_arrs[nbx](i, j, k) +
                        (dt * 0.5_rt *
                         (ddt_o_arrs[nbx](i, j, k) + ddt_arrs[nbx](i, j, k)));
                    // Defer n+1/2 update to time-stepping algorithm
                    // rho_nph(i, j, k) =
                    //     0.5_rt * (rho(i, j, k) + rho_o(i, j, k));
                });
        }
        amrex::Gpu::streamSynchronize();
    }

    void improve_explicit_diff(const amrex::Real /*unused*/) {}

    // data members
    PDEFields& fields;
};
} // namespace amr_wind::pde

#endif /* DENSITY_OPS_H */
