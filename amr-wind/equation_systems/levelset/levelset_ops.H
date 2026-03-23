#ifndef LEVELSET_OPS_H
#define LEVELSET_OPS_H

#include "amr-wind/equation_systems/levelset/levelset.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/** Field registration operator for Levelset equation system
 *  \ingroup levelset
 */
template <typename Scheme>
struct FieldRegOp<Levelset, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<Levelset, Scheme>(time, repo);

        auto& normal =
            repo.declare_cc_field("interface_normal", AMREX_SPACEDIM, 1, 1);
        auto& curvature = repo.declare_cc_field("interface_curvature", 1, 1, 1);

        normal.template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            repo.mesh(), time);
        curvature
            .template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
                repo.mesh(), time);

        // Register fields for output/restart
        auto& iomgr = sim.io_manager();
        iomgr.register_io_var(fields.field.name());
        iomgr.register_io_var(normal.name());
        iomgr.register_io_var(curvature.name());

        return fields;
    }

    CFDSim& sim;
};

/** Right-hand side (RHS) evaluator for Levelset transport equation
 *  \ingroup levelset
 */
template <typename Scheme>
struct ComputeRHSOp<Levelset, Scheme>
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
            const auto& phi_arrs = field(lev).arrays();
            const auto& phi_o_arrs = field_old(lev).const_arrays();
            const auto& ddt_o_arrs = conv_term(lev).const_arrays();

            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    phi_arrs[nbx](i, j, k) = phi_o_arrs[nbx](i, j, k) +
                                             (dt * ddt_o_arrs[nbx](i, j, k));
                    // Defer n+1/2 update to the time-stepping algorithm
                    // phi_nph(i, j, k) = 0.5_rt * (phi(i, j, k) + phi_o(i, j,
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
            auto phi_arrs = field(lev).arrays();
            const auto phi_o_arrs = field_old(lev).const_arrays();
            const auto ddt_arrs = conv_term(lev).const_arrays();
            const auto ddt_o_arrs = conv_term_old(lev).const_arrays();

            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    phi_arrs[nbx](i, j, k) =
                        phi_o_arrs[nbx](i, j, k) +
                        (dt * 0.5_rt *
                         (ddt_o_arrs[nbx](i, j, k) + ddt_arrs[nbx](i, j, k)));
                    // Defer n+1/2 update to time-stepping algorithm
                    // phi_nph(i, j, k) =
                    //     0.5_rt * (phi(i, j, k) + phi_o(i, j, k));
                });
        }
        amrex::Gpu::streamSynchronize();
    }

    void improve_explicit_diff(const amrex::Real /*unused*/) {}

    // data members
    PDEFields& fields;
};

/** Implementation of post-processing actions for the levelset equation
 *  \ingroup levelset
 */
template <>
struct PostSolveOp<Levelset>
{
    PostSolveOp(CFDSim& sim, PDEFields& fields) : m_sim(sim), m_fields(fields)
    {}

    void operator()(const amrex::Real /*unused*/)
    {
        m_fields.field.fillpatch(m_sim.time().current_time());
        auto& multiphase = m_sim.physics_manager().get<MultiPhase>();
        multiphase.set_density_via_levelset();
    }

    CFDSim& m_sim;
    PDEFields& m_fields;
};

} // namespace amr_wind::pde

#endif /* LEVELSET_OPS_H */
