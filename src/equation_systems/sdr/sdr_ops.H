#ifndef SDR_OPS_H
#define SDR_OPS_H

#include "amr-wind/equation_systems/sdr/SDR.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/** Implementation of Turbulence operator for SDR equation
 *  \ingroup sdr_eqn
 */
template <>
struct TurbulenceOp<SDR>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& mueff = m_fields.mueff;
        m_tmodel.update_scalar_diff(mueff, SDR::var_name());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

/** Implementation of post-processing actions for SDR equation
 *  \ingroup sdr_eqn
 */
template <>
struct PostSolveOp<SDR>
{
    PostSolveOp(CFDSim& /*unused*/, PDEFields& fields) : m_fields(fields) {}

    void operator()(const amrex::Real time)
    {
        field_ops::lower_bound(m_fields.field, clip_value);
        m_fields.field.fillpatch(time);
    }

    PDEFields& m_fields;
    amrex::Real clip_value{std::numeric_limits<float>::epsilon()};
};

template <typename Scheme>
struct FieldRegOp<SDR, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<SDR, Scheme>(time, repo);

        repo.declare_cc_field(
            SDR::var_name() + "_lhs_src_term", SDR::ndim, 1, 1);

        sim.io_manager().register_io_var(fields.field.name());

        return fields;
    }

    CFDSim& sim;
};

/** Diffusion operator for SDR equation
 *  \ingroup sdr_eqn
 */
template <typename Scheme>
struct DiffusionOp<SDR, Scheme> : public DiffSolverIface<typename SDR::MLDiffOp>
{
    static_assert(
        SDR::ndim == 1, "DiffusionOp invoked for non-scalar SDR type");
    static_assert(
        std::is_same_v<typename SDR::MLDiffOp, amrex::MLABecLaplacian>,
        "Invalid linear operator for scalar diffusion operator");

    DiffusionOp(
        PDEFields& fields, const bool has_overset, const bool mesh_mapping)
        : DiffSolverIface<typename SDR::MLDiffOp>(
              fields, has_overset, mesh_mapping)
        , m_lhs_src_term(
              fields.repo.get_field(SDR::var_name() + "_lhs_src_term"))
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::high));

        m_lhs_src_term.setVal(0.0_rt);
    }

    //! Computes the diffusion term that goes in the RHS
    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0_rt, -1.0_rt, fstate);

        auto tau_state =
            std::is_same_v<Scheme, fvm::Godunov> ? FieldState::New : fstate;
        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(
            this->m_pdefields.diff_term.state(tau_state).vec_ptrs(),
            this->m_pdefields.field.vec_ptrs());
    }

    void
    set_acoeffs(typename SDR::MLDiffOp& linop, const FieldState fstate) override
    {
        BL_PROFILE("amr-wind::pde::set_acoeffs");
        auto& repo = m_pdefields.repo;

        const int nlevels = repo.num_active_levels();
        auto& density = m_density.state(fstate);
        auto lhs_total = repo.create_scratch_field(1, 1);
        field_ops::lincomb(
            *lhs_total, 1.0_rt, m_lhs_src_term, 0, 1.0_rt, density, 0, 0, 1, 1);

        for (int lev = 0; lev < nlevels; ++lev) {
            linop.setACoeffs(lev, (*lhs_total)(lev));
        }
    }

    Field& m_lhs_src_term;
};

} // namespace amr_wind::pde

#endif /* SDR_OPS_H */
