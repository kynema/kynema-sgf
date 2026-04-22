#ifndef ICNS_OPS_H
#define ICNS_OPS_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/DiffusionOps.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

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
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(src_term, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& vf = src_term.array(mfi);
                const auto& rho = density(lev).const_array(mfi);
                const auto& gp = grad_p(lev).const_array(mfi);
                amrex::Array4<amrex::Real const> fac =
                    mesh_mapping ? ((*mesh_fac)(lev).const_array(mfi))
                                 : amrex::Array4<amrex::Real const>();

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        amrex::Real rhoinv = 1.0_rt / rho(i, j, k);
                        amrex::Real fac_x =
                            mesh_mapping ? (fac(i, j, k, 0)) : 1.0_rt;
                        amrex::Real fac_y =
                            mesh_mapping ? (fac(i, j, k, 1)) : 1.0_rt;
                        amrex::Real fac_z =
                            mesh_mapping ? (fac(i, j, k, 2)) : 1.0_rt;

                        vf(i, j, k, 0) =
                            -(1.0_rt / fac_x * gp(i, j, k, 0)) * rhoinv;
                        vf(i, j, k, 1) =
                            -(1.0_rt / fac_y * gp(i, j, k, 1)) * rhoinv;
                        vf(i, j, k, 2) =
                            -(1.0_rt / fac_z * gp(i, j, k, 2)) * rhoinv;
                    });

                for (const auto& src : this->sources) {
                    (*src)(lev, mfi, bx, fstate, vf);
                }
            }
        }
        // Multiply velocity src terms by density for momentum equation
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

} // namespace amr_wind::pde

#endif /* ICNS_OPS_H */
