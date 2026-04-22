#ifndef VOF_OPS_H
#define VOF_OPS_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind::pde {

/** Field registration operator for VOF equation system
 *  \ingroup vof
 */
template <typename Scheme>
struct FieldRegOp<VOF, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<VOF, Scheme>(time, repo);

        auto& levelset = repo.declare_cc_field("levelset", 1, 3, 1);
        auto& curvature = repo.declare_cc_field("interface_curvature", 1, 1, 1);

        levelset.set_default_fillpatch_bc(sim.time());
        curvature.set_default_fillpatch_bc(sim.time());

        // Register fields for output/restart
        auto& iomgr = sim.io_manager();
        iomgr.register_io_var(fields.field.name());

        return fields;
    }

    CFDSim& sim;
};

/** Right-hand side (RHS) evaluator for VOF transport equation
 */
template <typename Scheme>
struct ComputeRHSOp<VOF, Scheme>
{
    explicit ComputeRHSOp(PDEFields& fields_in) : fields(fields_in) {}

    void predictor_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real /*unused*/,
        bool /*unused*/)
    {}

    void corrector_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real /*unused*/,
        bool /*unused*/)
    {}

    void improve_explicit_diff(const amrex::Real /*unused*/) {}

    // data members
    PDEFields& fields;
};

/** Implementation of post-processing actions for the VOF equation
 *  \ingroup vof
 */

template <>
struct PostSolveOp<VOF>
{
    PostSolveOp(CFDSim& sim, PDEFields& fields) : m_sim(sim), m_fields(fields)
    {}

    void operator()(const amrex::Real /*unused*/)
    {
        m_fields.field.fillpatch(m_sim.time().new_time());
        auto& multiphase = m_sim.physics_manager().get<MultiPhase>();
        multiphase.set_density_via_vof();

        // Calculate n+1/2 density
        multiphase.set_nph_density();
        // Convert advalpha to advrho
        multiphase.calculate_advected_facedensity();
    }

    CFDSim& m_sim;
    PDEFields& m_fields;
};

} // namespace amr_wind::pde

#endif /* VOF_OPS_H */
