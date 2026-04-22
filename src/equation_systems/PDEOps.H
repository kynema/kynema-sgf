#ifndef PDEOPS_H
#define PDEOPS_H

#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/**
 *  \defgroup pdeop PDE Operators
 *  PDE Operators
 *
 *  The different terms in a PDE as well as other actions on a PDE are
 *  implemented within AMR-Wind using operator classes. These operators can be
 *  specialized for a given PDE or a Scheme based on specific traits that are
 *  passed to the amr_wind::pde::PDESytem class.
 *
 *  \ingroup eqsys
 */

/** %Field registration operator
 *  \ingroup pdeop
 *
 *  When a PDESystem is created, this operator is called first and is
 *  responsible for creating all the necessary fields that will be accessed by
 *  the PDE during the solution process.
 */
template <typename PDE, typename Scheme>
struct FieldRegOp
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    /** Perform initialization activities
     *
     *  Declares the necessary fields, sets their I/O properties
     *
     *  \return The collection of fields for this PDE
     */
    PDEFields operator()(const SimTime& time)
    {
        auto fields = create_fields_instance<PDE, Scheme>(time, sim.repo());

        // Register solution variable for this PDE as output/restart variable
        sim.io_manager().register_io_var(fields.field.name());
        return fields;
    }

    CFDSim& sim;
};

/** Base implementation of source terms for a PDE
 *  \ingroup pdeop
 *
 *  Holds a collection of source terms that are initialized using at runtime
 *  based on user inputs. The implementation consults the `source_terms` input
 *  variable within the PDE namespace (e.g., `ICNS.source_terms` for momentum
 *  sources) and initializes them. These source terms are then accumulated into
 *  the PDEFields::src_term during solution process.
 *
 *  \sa MomentumSource, TemperatureSource
 */
template <typename PDE>
struct SrcTermOpBase
{
    explicit SrcTermOpBase(PDEFields& fields_in)
        : fields(fields_in), m_density(fields_in.repo.get_field("density"))
    {}

    /** Initialize the source terms
     *
     *  At this stage the mesh has been created
     */
    void init_source_terms(const CFDSim& sim)
    {
        amrex::ParmParse pp(PDE::pde_name());
        amrex::Vector<std::string> src_terms;
        pp.queryarr("source_terms", src_terms);

        for (auto& src_name : src_terms) {
            // Prefer to use emplace_back here
            sources.emplace_back(PDE::SrcTerm::create(src_name, sim));
        }
    }

    //! Helper method to multiply the source terms with density
    void multiply_rho(const FieldState fstate)
    {
        BL_PROFILE("amr-wind::" + PDE::pde_name() + "::multiply_rho");
        const auto rhostate = field_impl::phi_state(fstate);
        const auto& density = m_density.state(rhostate);

        const int nlevels = fields.repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            auto& src_term = fields.src_term(lev);
            const auto& vf_arrs = src_term.arrays();
            const auto& rho_arrs = density(lev).const_arrays();

            amrex::ParallelFor(
                src_term, amrex::IntVect(0), fields.src_term.num_comp(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                    vf_arrs[nbx](i, j, k, n) *= rho_arrs[nbx](i, j, k);
                });
        }
        amrex::Gpu::streamSynchronize();
    }

    //! Update source terms during time-integration procedure
    virtual void
    operator()(const FieldState fstate, const bool /* mesh_mapping */)
    {
        // Zero out source term
        this->fields.src_term.setVal(0.0_rt);

        // Return early if there are no source terms to process
        if (this->sources.empty()) {
            return;
        }

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

                for (const auto& src : this->sources) {
                    (*src)(lev, mfi, bx, fstate, vf);
                }
            }
        }

        if (PDE::multiply_rho) {
            this->multiply_rho(fstate);
        }
    }

    PDEFields& fields;
    Field& m_density;
    amrex::Vector<std::unique_ptr<typename PDE::SrcTerm>> sources;
};

/** Implementation of source terms for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct SrcTermOp : SrcTermOpBase<PDE>
{
    explicit SrcTermOp(PDEFields& fields_in) : SrcTermOpBase<PDE>(fields_in) {}
};

template <typename PDE, typename Scheme, typename = void>
struct AdvectionOp
{};

template <typename PDE, typename Scheme, typename = void>
struct DiffusionOp
{};

/** Turbulence update operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct TurbulenceOp
{
    // cppcheck-suppress uninitMemberVar
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        m_tmodel.update_scalar_diff(m_fields.mueff, m_fields.field.name());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

/** Boundary condition update operator
 *  \ingroup pdeop
 */
template <typename PDE, typename = void>
struct BCOp;

/** Actions to be performed after a PDE solve
 *  \ingroup pdeop
 */
template <typename PDE>
struct PostSolveOp
{
    explicit PostSolveOp(CFDSim& sim, PDEFields& fields)
        : m_sim(sim), m_fields(fields)
    {}

    void operator()(const amrex::Real time) { m_fields.field.fillpatch(time); }

    CFDSim& m_sim;
    PDEFields& m_fields;
};

} // namespace amr_wind::pde

#endif /* PDEOPS_H */
