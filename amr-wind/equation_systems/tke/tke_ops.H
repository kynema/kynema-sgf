#ifndef TKE_OPS_H
#define TKE_OPS_H

#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/** Implementation of Turbulence operator for TKE equation
 *  \ingroup tke_eqn
 */
template <>
struct TurbulenceOp<TKE>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        auto& mueff = m_fields.mueff;
        m_tmodel.update_scalar_diff(mueff, TKE::var_name());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

/** Implementation of post-processing actions for TKE equation
 *  \ingroup tke_eqn
 */
template <>
struct PostSolveOp<TKE>
{
    PostSolveOp(CFDSim& /*unused*/, PDEFields& fields) : m_fields(fields) {}

    void operator()(const amrex::Real time)
    {
        // cppcheck-suppress constVariableReference
        auto& field = m_fields.field;
        const auto& repo = field.repo();
        const int nlevels = repo.num_active_levels();
        const auto clip_value = m_clip_value;
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& field_arrs = field(lev).arrays();
            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    field_arrs[nbx](i, j, k) =
                        (field_arrs[nbx](i, j, k) < 0.0_rt)
                            ? clip_value
                            : field_arrs[nbx](i, j, k);
                });
            const bool has_terrain = repo.int_field_exists("terrain_blank");
            if (has_terrain) {
                const IntField& m_terrain_blank =
                    repo.get_int_field("terrain_blank");
                const auto& blank_arrs = m_terrain_blank(lev).arrays();
                amrex::ParallelFor(
                    field(lev),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        field_arrs[nbx](i, j, k) =
                            (blank_arrs[nbx](i, j, k) == 1)
                                ? clip_value
                                : field_arrs[nbx](i, j, k);
                    });
            }
        }
        amrex::Gpu::streamSynchronize();
        m_fields.field.fillpatch(time);
    }

    PDEFields& m_fields;
    amrex::Real m_clip_value{std::numeric_limits<float>::epsilon()};
};

template <typename Scheme>
struct FieldRegOp<TKE, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in)
    {
        amrex::ParmParse pp("TKE");

        std::string interpolation;
        pp.query("interpolation", interpolation);

        if (interpolation == "PiecewiseConstant") {
            m_itype = FieldInterpolator::PiecewiseConstant;
            amrex::Print()
                << "using PiecewiseConstant interpolation for TKE fillpatch"
                << '\n';
        } else {
            amrex::Print() << "using default CellConsLinear interpolation for "
                              "TKE fillpatch"
                           << '\n';
        }
    }

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<TKE, Scheme>(time, repo, m_itype);

        repo.declare_cc_field(
            TKE::var_name() + "_lhs_src_term", TKE::ndim, 1, 1);

        sim.io_manager().register_io_var(fields.field.name());

        return fields;
    }

    CFDSim& sim;
    FieldInterpolator m_itype{FieldInterpolator::CellConsLinear};
};

/** Diffusion operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename Scheme>
struct DiffusionOp<TKE, Scheme> : public DiffSolverIface<typename TKE::MLDiffOp>
{
    static_assert(
        TKE::ndim == 1, "DiffusionOp invoked for non-scalar TKE type");
    static_assert(
        std::is_same_v<typename TKE::MLDiffOp, amrex::MLABecLaplacian>,
        "Invalid linear operator for scalar diffusion operator");

    DiffusionOp(
        PDEFields& fields, const bool has_overset, const bool mesh_mapping)
        : DiffSolverIface<typename TKE::MLDiffOp>(
              fields, has_overset, mesh_mapping)
        , m_lhs_src_term(
              fields.repo.get_field(TKE::var_name() + "_lhs_src_term"))
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
    set_acoeffs(typename TKE::MLDiffOp& linop, const FieldState fstate) override
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

#endif /* TKE_OPS_H */
