#ifndef ICNS_OPS_H
#define ICNS_OPS_H

#include "src/equation_systems/PDETraits.H"
#include "src/equation_systems/PDEOps.H"
#include "src/equation_systems/PDEHelpers.H"
#include "src/equation_systems/AdvOp_Godunov.H"
#include "src/equation_systems/AdvOp_MOL.H"
#include "src/equation_systems/DiffusionOps.H"
#include "src/equation_systems/icns/icns.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde {

/** Specialization of field registration for ICNS
 *  \ingroup icns
 *
 *  In addition to the PDEFields, this operator initializes several additional
 *  fields necessary to solve the Navier-Stokes equation system.
 *
 *  - `density` - The fluid density field variable
 *  - `pressure` - The perturbation pressure field
 *  - `umac, vmac, wmac` - The face velocity fields used for flux computations
 *  - `grad_p` - Pressure gradient used for projection
 */
template <typename Scheme>
struct FieldRegOp<ICNS, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<ICNS, Scheme>(time, repo);

        auto& rho = repo.declare_cc_field(
            "density", 1, Scheme::nghost_state, Scheme::num_states);
        auto& grad_p = repo.declare_cc_field("gp", ICNS::ndim, 0, 1);
        auto& pressure = repo.declare_nd_field("p", 1, Scheme::nghost_state, 1);
        repo.declare_face_normal_field(
            {"u_mac", "v_mac", "w_mac"}, 1, Scheme::nghost_mac, 1);

        rho.template register_fill_patch_op<
            FieldFillPatchOps<FieldBCDirichlet>>(repo.mesh(), time);
        grad_p.template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            repo.mesh(), time);
        // pressure.template
        // register_fill_patch_op<FieldFillConstScalar>(0.0_rt); fixme this is
        // only necessary since tioga does not fill in ghosts convert back later
        pressure
            .template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
                repo.mesh(), time);
        rho.fillpatch_on_regrid() = true;
        grad_p.fillpatch_on_regrid() = true;

        // Density always requires half state
        rho.create_state(FieldState::NPH);

        // Register fields for output/restart
        auto& iomgr = sim.io_manager();
        iomgr.register_io_var(fields.field.name());
        iomgr.register_io_var(rho.name());
        iomgr.register_io_var(pressure.name());
        iomgr.register_io_var(grad_p.name());
        return fields;
    }

    CFDSim& sim;
};

/** Specialization of the source term operator for ICNS
 *  \ingroup icns
 */
template <>
struct SrcTermOp<ICNS> : SrcTermOpBase<ICNS>
{
    explicit SrcTermOp(PDEFields& fields_in)
        : SrcTermOpBase<ICNS>(fields_in), grad_p(fields_in.repo.get_field("gp"))
    {}

    void operator()(const FieldState fstate, const bool mesh_mapping) override
    {
        const auto rhostate = field_impl::phi_state(fstate);
        const auto& density = m_density.state(rhostate);
        Field const* mesh_fac =
            mesh_mapping
                ? &(this->fields.repo.get_mesh_mapping_field(FieldLoc::CELL))
                : nullptr;

        const int nlevels = this->fields.repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            auto& src_term = this->fields.src_term(lev);

            auto const& src_arrs = src_term.arrays();
            auto const& rho_arrs = density(lev).const_arrays();
            auto const& gp_arrs = grad_p(lev).const_arrays();
            auto const& fac_arrs =
                mesh_mapping ? (*mesh_fac)(lev).const_arrays()
                             : amrex::MultiArray4<amrex::Real const>();

            amrex::ParallelFor(
                src_term, amrex::IntVect(0), AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                    const amrex::Real rhoinv = 1.0_rt / rho_arrs[nbx](i, j, k);
                    const amrex::Real fac_n =
                        mesh_mapping ? fac_arrs[nbx](i, j, k, n) : 1.0_rt;
                    src_arrs[nbx](i, j, k, n) =
                        -(1.0_rt / fac_n * gp_arrs[nbx](i, j, k, n)) * rhoinv;
                });
            amrex::Gpu::streamSynchronize();

            for (const auto& src : this->sources) {
                (*src)(lev, fstate, src_term);
            }
            amrex::Gpu::streamSynchronize();
        }
        this->multiply_rho(fstate);
    }

    Field& grad_p;
};

/** Effective turbulent viscosity computation for ICNS
 *  \ingroup icns
 */
template <>
struct TurbulenceOp<ICNS>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()() { m_tmodel.update_mueff(m_fields.mueff); }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

} // namespace kynema_sgf::pde

#endif /* ICNS_OPS_H */
