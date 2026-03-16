#include "amr-wind/equation_systems/DiffusionOps.H"
#include "amr-wind/utilities/console_io.H"

#include "AMReX_MLTensorOp.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

template <typename LinOp>
DiffSolverIface<LinOp>::DiffSolverIface(
    PDEFields& fields,
    const bool has_overset,
    const bool mesh_mapping,
    const std::string& prefix)
    : m_pdefields(fields)
    , m_density(fields.repo.get_field("density"))
    , m_options(prefix, m_pdefields.field.name() + "_" + prefix)
    , m_mesh_mapping(mesh_mapping)
{
    amrex::LPInfo isolve = m_options.lpinfo();
    amrex::LPInfo iapply;

    iapply.setMaxCoarseningLevel(0);

    const auto& mesh = m_pdefields.repo.mesh();
    if (!has_overset) {
        m_solver.reset(new LinOp(
            mesh.Geom(0, mesh.finestLevel()),
            mesh.boxArray(0, mesh.finestLevel()),
            mesh.DistributionMap(0, mesh.finestLevel()), isolve));
        m_applier.reset(new LinOp(
            mesh.Geom(0, mesh.finestLevel()),
            mesh.boxArray(0, mesh.finestLevel()),
            mesh.DistributionMap(0, mesh.finestLevel()), iapply));
    } else {
        auto imask = fields.repo.get_int_field("mask_cell").vec_const_ptrs();
        m_solver.reset(new LinOp(
            mesh.Geom(0, mesh.finestLevel()),
            mesh.boxArray(0, mesh.finestLevel()),
            mesh.DistributionMap(0, mesh.finestLevel()), imask, isolve));
        m_applier.reset(new LinOp(
            mesh.Geom(0, mesh.finestLevel()),
            mesh.boxArray(0, mesh.finestLevel()),
            mesh.DistributionMap(0, mesh.finestLevel()), imask, iapply));
    }

    m_solver->setMaxOrder(m_options.max_order);
    m_applier->setMaxOrder(m_options.max_order);

    // It is the sub-classes responsibility to set the linear solver BC for the
    // operators.
}

template <typename LinOp>
void DiffSolverIface<LinOp>::setup_operator(
    LinOp& linop,
    const amrex::Real alpha,
    const amrex::Real beta,
    const FieldState fstate)
{
    BL_PROFILE("amr-wind::setup_operator");
    auto& repo = m_pdefields.repo;
    const int nlevels = repo.num_active_levels();

    linop.setScalars(alpha, beta);
    for (int lev = 0; lev < nlevels; ++lev) {
        linop.setLevelBC(lev, &m_pdefields.field(lev));
    }
    this->set_acoeffs(linop, fstate);
    set_bcoeffs(linop);
}

template <typename LinOp>
void DiffSolverIface<LinOp>::set_acoeffs(LinOp& linop, const FieldState fstate)
{
    BL_PROFILE("amr-wind::set_acoeffs");
    auto& repo = m_pdefields.repo;
    const int nlevels = repo.num_active_levels();
    const auto& density = m_density.state(fstate);
    Field const* mesh_detJ =
        m_mesh_mapping ? &(repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                       : nullptr;
    std::unique_ptr<ScratchField> rho_times_detJ =
        m_mesh_mapping ? repo.create_scratch_field(
                             1, m_density.num_grow()[0], FieldLoc::CELL)
                       : nullptr;

    for (int lev = 0; lev < nlevels; ++lev) {
        if (m_mesh_mapping) {
            (*rho_times_detJ)(lev).setVal(0.0_rt);
            amrex::MultiFab::AddProduct(
                (*rho_times_detJ)(lev), density(lev), 0, (*mesh_detJ)(lev), 0,
                0, 1, m_density.num_grow()[0]);
            linop.setACoeffs(lev, (*rho_times_detJ)(lev));
        } else {
            linop.setACoeffs(lev, density(lev));
        }
    }
}

template <typename LinOp>
void DiffSolverIface<LinOp>::setup_solver(amrex::MLMG& mlmg)
{
    BL_PROFILE("amr-wind::setup_solver");
    // Set all MLMG options
    m_options(mlmg);
}

template <typename LinOp>
void DiffSolverIface<LinOp>::linsys_solve_impl()
{
    FieldState fstate = FieldState::New;
    auto& repo = this->m_pdefields.repo;
    auto& field = this->m_pdefields.field;
    if (field.in_uniform_space()) {
        amrex::Abort(
            "For diffusion solve, velocity should not be in uniform mesh "
            "space.");
    }
    const auto& density = m_density.state(fstate);
    const int nlevels = repo.num_active_levels();
    const int ndim = field.num_comp();
    auto rhs_ptr = repo.create_scratch_field("rhs", field.num_comp(), 0);

    // Always multiply with rho since there is no diffusion term for density
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& rhs = (*rhs_ptr)(lev);
        const auto& rhs_arrs = rhs.arrays();
        const auto& fld_arrs = field(lev).const_arrays();
        const auto& rho_arrs = density(lev).const_arrays();

        amrex::ParallelFor(
            rhs, amrex::IntVect(0), ndim,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                rhs_arrs[nbx](i, j, k, n) =
                    rho_arrs[nbx](i, j, k) * fld_arrs[nbx](i, j, k, n);
            });
    }
    amrex::Gpu::streamSynchronize();

    amrex::MLMG mlmg(*this->m_solver);
    this->setup_solver(mlmg);

    mlmg.solve(
        field.vec_ptrs(), rhs_ptr->vec_const_ptrs(), this->m_options.rel_tol,
        this->m_options.abs_tol);

    io::print_mlmg_info(field.name() + "_solve", mlmg);
}

template <typename LinOp>
void DiffSolverIface<LinOp>::linsys_solve(const amrex::Real dt)
{
    FieldState fstate = FieldState::New;
    this->setup_operator(*this->m_solver, 1.0_rt, dt, fstate);
    this->linsys_solve_impl();
}

template class DiffSolverIface<amrex::MLABecLaplacian>;
template class DiffSolverIface<amrex::MLTensorOp>;

} // namespace amr_wind::pde
