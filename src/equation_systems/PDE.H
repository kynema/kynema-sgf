#ifndef PDE_H
#define PDE_H

#include <string>
#include "amr-wind/CFDSim.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/CompRHSOps.H"
#include "amr-wind/equation_systems/DiffusionOps.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/** Implementation of the PDE interface for transport equations
 *  \ingroup eqsys
 *
 *  This class takes two traits, one describing the PDE type and the other
 *  describing the numerical scheme (e.g., Godunov, MOL) to generate a family of
 *  PDE systems and its solution process. Each operation is performed by a
 *  separate operator object that can be customized based on the different
 *  traits as necessary.
 *
 *  \tparam PDE The PDE trait (e.g., ICNS, Temperature, etc.)
 *  \tparam Scheme Numerical scheme trait (e.g., Godunov, MOL)
 */
template <typename PDE, typename Scheme>
class PDESystem : public PDEBase::Register<PDESystem<PDE, Scheme>>
{
public:
    using PDEType = PDE;
    using SchemeType = Scheme;

    //! Unique identifier used to register and create this instance on-demand
    //! through the factor interface
    static std::string identifier()
    {
        return PDE::pde_name() + "-" + Scheme::scheme_name();
    }

    /**
     *  \param sim A CFDSim instance
     */
    explicit PDESystem(CFDSim& sim)
        : m_sim(sim)
        , m_time(sim.time())
        , m_repo(sim.repo())
        , m_fields(FieldRegOp<PDE, Scheme>(m_sim)(m_time))
        , m_src_op(m_fields)
        , m_rhs_op(m_fields)
        , m_bc_op(m_fields, m_time)
        , m_post_solve_op(m_sim, m_fields)
    {
        m_bc_op.init_bcs();
    }

    //! Perform initialization actions after the mesh is created
    void initialize() override
    {
        if (PDE::has_diffusion) {
            BL_PROFILE("amr-wind::" + this->identifier() + "::initialize");
            m_diff_op.reset(new DiffusionOp<PDE, Scheme>(
                m_fields, m_sim.has_overset(), m_sim.has_mesh_mapping()));
            m_turb_op.reset(
                new TurbulenceOp<PDE>(m_sim.turbulence_model(), m_fields));
        }

        const bool variable_density =
            (!m_sim.pde_manager().constant_density() ||
             m_sim.physics_manager().contains("MultiPhase"));

        m_adv_op.reset(new AdvectionOp<PDE, Scheme>(
            m_sim, m_fields, m_sim.has_overset(), variable_density,
            m_sim.has_mesh_mapping(), m_sim.is_anelastic()));
        m_src_op.init_source_terms(m_sim);

        // Post-solve operations should also apply after initialization
        m_post_solve_op(m_time.current_time());
    }

    //! Perform update actions after a regrid is performed
    void post_regrid_actions() override
    {
        if (PDE::has_diffusion) {
            BL_PROFILE(
                "amr-wind::" + this->identifier() + "::post_regrid_actions");
            m_diff_op.reset(new DiffusionOp<PDE, Scheme>(
                m_fields, m_sim.has_overset(), m_sim.has_mesh_mapping()));
        }

        const bool variable_density =
            (!m_sim.pde_manager().constant_density() ||
             m_sim.physics_manager().contains("MultiPhase"));

        m_adv_op.reset(new AdvectionOp<PDE, Scheme>(
            m_sim, m_fields, m_sim.has_overset(), variable_density,
            m_sim.has_mesh_mapping(), m_sim.is_anelastic()));

        // Post-solve operations should also apply after a regrid
        m_post_solve_op(m_time.current_time());
    }

    //! Return the object holding the fields necessary for solving this PDE
    PDEFields& fields() override { return m_fields; }
    [[nodiscard]] const PDEFields& fields() const override { return m_fields; }

    void compute_source_term(const FieldState fstate) override
    {
        BL_PROFILE("amr-wind::" + this->identifier() + "::compute_source_term");
        m_src_op(fstate, m_sim.has_mesh_mapping());
    }

    void compute_mueff(const FieldState /*unused*/) override
    {
        if (PDE::has_diffusion) {
            BL_PROFILE("amr-wind::" + this->identifier() + "::compute_mueff");
            (*m_turb_op)();
        }
    }

    void compute_diffusion_term(const FieldState fstate) override
    {
        if (PDE::has_diffusion) {
            BL_PROFILE(
                "amr-wind::" + this->identifier() + "::compute_diffusion_term");
            m_bc_op.apply_bcs(fstate);
            m_diff_op->compute_diff_term(fstate);
        }
    }

    void compute_advection_term(const FieldState fstate) override
    {
        BL_PROFILE(
            "amr-wind::" + this->identifier() + "::compute_advection_term");
        (*m_adv_op)(fstate, m_time.delta_t());
    }

    void pre_advection_actions(const FieldState fstate) override
    {
        BL_PROFILE(
            "amr-wind::" + this->identifier() + "::pre_advection_actions");
        m_adv_op->preadvect(
            fstate, m_time.delta_t(),
            0.5_rt * (m_time.current_time() + m_time.new_time()));
    }

    void compute_predictor_rhs(const DiffusionType difftype) override
    {
        BL_PROFILE(
            "amr-wind::" + this->identifier() + "::compute_predictor_rhs");
        m_rhs_op.predictor_rhs(
            difftype, m_time.delta_t(), m_sim.has_mesh_mapping());
    }

    void compute_corrector_rhs(const DiffusionType difftype) override
    {
        BL_PROFILE(
            "amr-wind::" + this->identifier() + "::compute_corrector_rhs");
        m_rhs_op.corrector_rhs(
            difftype, m_time.delta_t(), m_sim.has_mesh_mapping());
    }

    void solve(const amrex::Real dt) override
    {
        if (PDE::has_diffusion) {
            BL_PROFILE("amr-wind::" + this->identifier() + "::linsys_solve");
            m_bc_op.apply_bcs(FieldState::New);
            m_diff_op->linsys_solve(dt);
        }
    }

    void post_solve_actions() override { m_post_solve_op(m_time.new_time()); }

    void improve_explicit_diffusion(const amrex::Real dt) override
    {
        if (PDE::has_diffusion) {
            BL_PROFILE(
                "amr-wind::" + this->identifier() +
                "::improve_explicit_diffusion");
            m_rhs_op.improve_explicit_diff(dt);
        }
    }

protected:
    //! CFD simulation controller instance
    CFDSim& m_sim;

    //! Time controls instance
    const SimTime& m_time;

    //! Field repository
    FieldRepo& m_repo;

    PDEFields m_fields;

    //! Source term computation operator
    SrcTermOp<PDE> m_src_op;

    //! Advection term computation operator
    std::unique_ptr<AdvectionOp<PDE, Scheme>> m_adv_op;

    //! RHS computation operator
    ComputeRHSOp<PDE, Scheme> m_rhs_op;

    //! Diffusion term computation operator
    std::unique_ptr<DiffusionOp<PDE, Scheme>> m_diff_op;

    //! Turbulence operator
    std::unique_ptr<TurbulenceOp<PDE>> m_turb_op;

    BCOp<PDE> m_bc_op;

    PostSolveOp<PDE> m_post_solve_op;
};

} // namespace amr_wind::pde

#endif /* PDE_H */
