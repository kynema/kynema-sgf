#ifndef ICNS_DIFFUSION_H
#define ICNS_DIFFUSION_H

#include <memory>

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/DiffusionOps.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/utilities/console_io.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

class ICNSDiffTensorOp : public DiffSolverIface<amrex::MLTensorOp>
{
public:
    ICNSDiffTensorOp(
        PDEFields& fields, const bool has_overset, const bool mesh_mapping)
        : DiffSolverIface<amrex::MLTensorOp>(fields, has_overset, mesh_mapping)
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::high));
    }

    template <typename Scheme>
    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0_rt, -1.0_rt, fstate);

        auto tau_state =
            std::is_same_v<Scheme, fvm::Godunov> ? FieldState::New : fstate;
        auto& divtau = this->m_pdefields.diff_term.state(tau_state);

        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(divtau.vec_ptrs(), this->m_pdefields.field.vec_ptrs());
    }
};

class ICNSDiffScalarOp
{
public:
    ICNSDiffScalarOp(
        PDEFields& fields,
        const bool has_overset,
        const bool mesh_mapping,
        const std::string& prefix = "diffusion")
        : m_pdefields(fields)
        , m_density(fields.repo.get_field("density"))
        , m_options(prefix, m_pdefields.field.name() + "_" + prefix)
        , m_mesh_mapping(mesh_mapping)
    {
        amrex::LPInfo isolve = m_options.lpinfo();
        amrex::LPInfo iapply;

        iapply.setMaxCoarseningLevel(0);
        const auto& mesh = m_pdefields.repo.mesh();

        const auto& bclo = diffusion::get_diffuse_tensor_bc(
            m_pdefields.field, amrex::Orientation::low);
        const auto& bchi = diffusion::get_diffuse_tensor_bc(
            m_pdefields.field, amrex::Orientation::high);

        if (!has_overset) {

            m_solver_scalar = std::make_unique<amrex::MLABecLaplacian>(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), isolve,
                amrex::Vector<
                    amrex::FabFactory<amrex::MLABecLaplacian::FAB> const*>(),
                AMREX_SPACEDIM);
            m_applier_scalar = std::make_unique<amrex::MLABecLaplacian>(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), iapply,
                amrex::Vector<
                    amrex::FabFactory<amrex::MLABecLaplacian::FAB> const*>(),
                AMREX_SPACEDIM);
        } else {
            auto imask =
                fields.repo.get_int_field("mask_cell").vec_const_ptrs();
            m_solver_scalar = std::make_unique<amrex::MLABecLaplacian>(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), imask, isolve,
                amrex::Vector<
                    amrex::FabFactory<amrex::MLABecLaplacian::FAB> const*>(),
                AMREX_SPACEDIM);
            m_applier_scalar = std::make_unique<amrex::MLABecLaplacian>(
                mesh.Geom(0, mesh.finestLevel()),
                mesh.boxArray(0, mesh.finestLevel()),
                mesh.DistributionMap(0, mesh.finestLevel()), imask, iapply,
                amrex::Vector<
                    amrex::FabFactory<amrex::MLABecLaplacian::FAB> const*>(),
                AMREX_SPACEDIM);
        }

        m_solver_scalar->setMaxOrder(m_options.max_order);
        m_applier_scalar->setMaxOrder(m_options.max_order);

        m_solver_scalar->setDomainBC(bclo, bchi);
        m_applier_scalar->setDomainBC(bclo, bchi);
    }

    template <typename Scheme>
    void compute_diff_term(const FieldState fstate)
    {
        auto tau_state =
            std::is_same_v<Scheme, fvm::Godunov> ? FieldState::New : fstate;
        const auto& repo = m_pdefields.repo;
        const int nlevels = repo.num_active_levels();
        const auto& geom = repo.mesh().Geom();

        bool diff_for_RHS(fstate == amr_wind::FieldState::New);
        auto& divtau = m_pdefields.diff_term.state(tau_state);
        const auto& density = m_density.state(fstate);
        const auto& viscosity = m_pdefields.mueff;
        Field const* mesh_detJ =
            m_mesh_mapping ? &(repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                           : nullptr;
        std::unique_ptr<ScratchField> rho_times_detJ =
            m_mesh_mapping ? repo.create_scratch_field(
                                 1, m_density.num_grow()[0], FieldLoc::CELL)
                           : nullptr;

        const amrex::Real alpha = 0.0_rt;
        const amrex::Real beta = -1.0_rt;
        m_applier_scalar->setScalars(alpha, beta);

        for (int lev = 0; lev < nlevels; ++lev) {
            m_applier_scalar->setLevelBC(lev, &m_pdefields.field(lev));

            // A coeffs
            if (m_mesh_mapping) {
                (*rho_times_detJ)(lev).setVal(0.0_rt);
                amrex::MultiFab::AddProduct(
                    (*rho_times_detJ)(lev), density(lev), 0, (*mesh_detJ)(lev),
                    0, 0, 1, m_density.num_grow()[0]);
                m_applier_scalar->setACoeffs(lev, (*rho_times_detJ)(lev));
            } else {
                m_applier_scalar->setACoeffs(lev, density(lev));
            }

            // B coeffs
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));
            if (m_mesh_mapping) {
                diffusion::viscosity_to_uniform_space(b, repo, lev);
            }
            m_applier_scalar->setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }

        amrex::MLMG mlmg(*m_applier_scalar);
        mlmg.apply(divtau.vec_ptrs(), m_pdefields.field.vec_ptrs());

        if (!diff_for_RHS) {
            for (int lev = 0; lev < nlevels; ++lev) {
                const auto& divtau_arrs = divtau(lev).arrays();
                const auto& rho_arrs = density(lev).const_arrays();

                amrex::ParallelFor(
                    divtau(lev),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        amrex::Real rhoinv = 1.0_rt / rho_arrs[nbx](i, j, k);
                        divtau_arrs[nbx](i, j, k, 0) *= rhoinv;
                        divtau_arrs[nbx](i, j, k, 1) *= rhoinv;
                        divtau_arrs[nbx](i, j, k, 2) *= rhoinv;
                    });
            }
            amrex::Gpu::streamSynchronize();
        }
    }

    void linsys_solve(const amrex::Real dt)
    {
        const FieldState fstate = FieldState::New;
        auto& repo = m_pdefields.repo;
        const auto& geom = repo.mesh().Geom();
        const auto& field = m_pdefields.field;
        const auto& density = m_density.state(fstate);
        const int nlevels = repo.num_active_levels();
        const int ndim = field.num_comp();
        auto rhs_ptr = repo.create_scratch_field("rhs", field.num_comp(), 0);
        const auto& viscosity = m_pdefields.mueff;
        Field const* mesh_detJ =
            m_mesh_mapping ? &(repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                           : nullptr;
        std::unique_ptr<ScratchField> rho_times_detJ =
            m_mesh_mapping ? repo.create_scratch_field(
                                 1, m_density.num_grow()[0], FieldLoc::CELL)
                           : nullptr;

        const amrex::Real alpha = 1.0_rt;
        const amrex::Real beta = dt;
        m_solver_scalar->setScalars(alpha, beta);

        for (int lev = 0; lev < nlevels; ++lev) {
            m_solver_scalar->setLevelBC(lev, &m_pdefields.field(lev));

            // A coeffs
            if (m_mesh_mapping) {
                (*rho_times_detJ)(lev).setVal(0.0_rt);
                amrex::MultiFab::AddProduct(
                    (*rho_times_detJ)(lev), density(lev), 0, (*mesh_detJ)(lev),
                    0, 0, 1, m_density.num_grow()[0]);
                m_solver_scalar->setACoeffs(lev, (*rho_times_detJ)(lev));
            } else {
                m_solver_scalar->setACoeffs(lev, density(lev));
            }

            // B coeffs
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));
            if (m_mesh_mapping) {
                diffusion::viscosity_to_uniform_space(b, repo, lev);
            }
            m_solver_scalar->setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }

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

        amrex::MLMG mlmg(*m_solver_scalar);
        m_options(mlmg);
        mlmg.solve(
            m_pdefields.field.vec_ptrs(), rhs_ptr->vec_const_ptrs(),
            m_options.rel_tol, m_options.abs_tol);

        io::print_mlmg_info(field.name() + "_multicomponent_solve", mlmg);
    }

protected:
    PDEFields& m_pdefields;
    Field& m_density;
    MLMGOptions m_options;
    bool m_mesh_mapping{false};

    std::unique_ptr<amrex::MLABecLaplacian> m_solver_scalar;
    std::unique_ptr<amrex::MLABecLaplacian> m_applier_scalar;
};

class ICNSDiffScalarSegregatedOp
{
public:
    ICNSDiffScalarSegregatedOp(
        PDEFields& fields,
        const bool has_overset,
        const bool mesh_mapping,
        const std::string& prefix = "diffusion")
        : m_pdefields(fields)
        , m_density(fields.repo.get_field("density"))
        , m_options(prefix, m_pdefields.field.name() + "_" + prefix)
        , m_mesh_mapping(mesh_mapping)
    {
        amrex::LPInfo isolve = m_options.lpinfo();
        amrex::LPInfo iapply;

        iapply.setMaxCoarseningLevel(0);
        const auto& mesh = m_pdefields.repo.mesh();

        const auto& bclo = diffusion::get_diffuse_tensor_bc(
            m_pdefields.field, amrex::Orientation::low);
        const auto& bchi = diffusion::get_diffuse_tensor_bc(
            m_pdefields.field, amrex::Orientation::high);

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (!has_overset) {
                m_solver_scalar[i] = std::make_unique<amrex::MLABecLaplacian>(
                    mesh.Geom(0, mesh.finestLevel()),
                    mesh.boxArray(0, mesh.finestLevel()),
                    mesh.DistributionMap(0, mesh.finestLevel()), isolve);
                m_applier_scalar[i] = std::make_unique<amrex::MLABecLaplacian>(
                    mesh.Geom(0, mesh.finestLevel()),
                    mesh.boxArray(0, mesh.finestLevel()),
                    mesh.DistributionMap(0, mesh.finestLevel()), iapply);
            } else {
                auto imask =
                    fields.repo.get_int_field("mask_cell").vec_const_ptrs();
                m_solver_scalar[i] = std::make_unique<amrex::MLABecLaplacian>(
                    mesh.Geom(0, mesh.finestLevel()),
                    mesh.boxArray(0, mesh.finestLevel()),
                    mesh.DistributionMap(0, mesh.finestLevel()), imask, isolve);
                m_applier_scalar[i] = std::make_unique<amrex::MLABecLaplacian>(
                    mesh.Geom(0, mesh.finestLevel()),
                    mesh.boxArray(0, mesh.finestLevel()),
                    mesh.DistributionMap(0, mesh.finestLevel()), imask, iapply);
            }

            m_solver_scalar[i]->setMaxOrder(m_options.max_order);
            m_applier_scalar[i]->setMaxOrder(m_options.max_order);

            m_solver_scalar[i]->setDomainBC(bclo[i], bchi[i]);
            m_applier_scalar[i]->setDomainBC(bclo[i], bchi[i]);
        }
    }

    template <typename Scheme>
    void compute_diff_term(const FieldState fstate)
    {
        auto tau_state =
            std::is_same_v<Scheme, fvm::Godunov> ? FieldState::New : fstate;
        const auto& repo = m_pdefields.repo;
        const int nlevels = repo.num_active_levels();
        const auto& geom = repo.mesh().Geom();

        auto& divtau = m_pdefields.diff_term.state(tau_state);
        bool diff_for_RHS(fstate == amr_wind::FieldState::New);
        const auto& density = m_density.state(fstate);
        const auto& viscosity = m_pdefields.mueff;
        Field const* mesh_detJ =
            m_mesh_mapping ? &(repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                           : nullptr;
        std::unique_ptr<ScratchField> rho_times_detJ =
            m_mesh_mapping ? repo.create_scratch_field(
                                 1, m_density.num_grow()[0], FieldLoc::CELL)
                           : nullptr;

        const amrex::Real alpha = 0.0_rt;
        const amrex::Real beta = -1.0_rt;

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            m_applier_scalar[i]->setScalars(alpha, beta);

            for (int lev = 0; lev < nlevels; ++lev) {

                m_applier_scalar[i]->setLevelBC(
                    lev, &m_pdefields.field.subview(i)(lev));

                if (m_mesh_mapping) {
                    (*rho_times_detJ)(lev).setVal(0.0_rt);
                    amrex::MultiFab::AddProduct(
                        (*rho_times_detJ)(lev), density(lev), 0,
                        (*mesh_detJ)(lev), 0, 0, 1, m_density.num_grow()[0]);
                }

                // A coeffs
                if (m_mesh_mapping) {
                    m_applier_scalar[i]->setACoeffs(
                        lev, (*rho_times_detJ)(lev));
                } else {
                    m_applier_scalar[i]->setACoeffs(lev, density(lev));
                }

                // B coeffs
                auto b = diffusion::average_velocity_eta_to_faces(
                    geom[lev], viscosity(lev));

                if (m_mesh_mapping) {
                    diffusion::viscosity_to_uniform_space(b, repo, lev);
                }

                m_applier_scalar[i]->setBCoeffs(
                    lev, amrex::GetArrOfConstPtrs(b));
            }

            auto divtau_comp = divtau.subview(i);
            auto vel_comp = m_pdefields.field.subview(i);

            amrex::MLMG mlmg(*m_applier_scalar[i]);
            mlmg.apply(divtau_comp.vec_ptrs(), vel_comp.vec_ptrs());
        }

        if (!diff_for_RHS) {
            for (int lev = 0; lev < nlevels; ++lev) {
                const auto& divtau_arrs = divtau(lev).arrays();
                const auto& rho_arrs = density(lev).const_arrays();

                amrex::ParallelFor(
                    divtau(lev),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        amrex::Real rhoinv = 1.0_rt / rho_arrs[nbx](i, j, k);
                        divtau_arrs[nbx](i, j, k, 0) *= rhoinv;
                        divtau_arrs[nbx](i, j, k, 1) *= rhoinv;
                        divtau_arrs[nbx](i, j, k, 2) *= rhoinv;
                    });
            }
            amrex::Gpu::streamSynchronize();
        }
    }

    void linsys_solve(const amrex::Real dt)
    {
        const FieldState fstate = FieldState::New;
        auto& repo = m_pdefields.repo;
        const auto& field = m_pdefields.field;
        const auto& density = m_density.state(fstate);
        const int nlevels = repo.num_active_levels();
        const int ndim = field.num_comp();
        auto rhs_ptr = repo.create_scratch_field("rhs", field.num_comp(), 0);
        const auto& viscosity = m_pdefields.mueff;
        const auto& geom = repo.mesh().Geom();

        Field const* mesh_detJ =
            m_mesh_mapping ? &(repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                           : nullptr;
        std::unique_ptr<ScratchField> rho_times_detJ =
            m_mesh_mapping ? repo.create_scratch_field(
                                 1, m_density.num_grow()[0], FieldLoc::CELL)
                           : nullptr;

        const amrex::Real alpha = 1.0_rt;
        const amrex::Real beta = dt;

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            m_solver_scalar[i]->setScalars(alpha, beta);

            for (int lev = 0; lev < nlevels; ++lev) {

                if (m_mesh_mapping) {
                    (*rho_times_detJ)(lev).setVal(0.0_rt);
                    amrex::MultiFab::AddProduct(
                        (*rho_times_detJ)(lev), density(lev), 0,
                        (*mesh_detJ)(lev), 0, 0, 1, m_density.num_grow()[0]);
                }

                m_solver_scalar[i]->setLevelBC(
                    lev, &m_pdefields.field.subview(i)(lev));

                // A coeffs
                if (m_mesh_mapping) {
                    m_solver_scalar[i]->setACoeffs(lev, (*rho_times_detJ)(lev));
                } else {
                    m_solver_scalar[i]->setACoeffs(lev, density(lev));
                }

                // B coeffs
                auto b = diffusion::average_velocity_eta_to_faces(
                    geom[lev], viscosity(lev));
                if (m_mesh_mapping) {
                    diffusion::viscosity_to_uniform_space(b, repo, lev);
                }

                m_solver_scalar[i]->setBCoeffs(
                    lev, amrex::GetArrOfConstPtrs(b));
            }
        }

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

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            auto vel_comp = m_pdefields.field.subview(i);
            auto rhs_ptr_comp = rhs_ptr->subview(i);

            amrex::MLMG mlmg(*m_solver_scalar[i]);
            m_options(mlmg);

            mlmg.solve(
                vel_comp.vec_ptrs(), rhs_ptr_comp.vec_const_ptrs(),
                m_options.rel_tol, m_options.abs_tol);

            io::print_mlmg_info(
                field.name() + std::to_string(i) + "_solve", mlmg);
        }
    }

protected:
    PDEFields& m_pdefields;
    Field& m_density;
    MLMGOptions m_options;
    bool m_mesh_mapping{false};

    amrex::Array<std::unique_ptr<amrex::MLABecLaplacian>, AMREX_SPACEDIM>
        m_solver_scalar;
    amrex::Array<std::unique_ptr<amrex::MLABecLaplacian>, AMREX_SPACEDIM>
        m_applier_scalar;
};

/** Specialization of diffusion operator for ICNS
 *  \ingroup icns
 */
template <typename Scheme>
struct DiffusionOp<ICNS, Scheme>
{
    std::unique_ptr<ICNSDiffTensorOp> m_tensor_op;
    std::unique_ptr<ICNSDiffScalarOp> m_scalar_op;
    std::unique_ptr<ICNSDiffScalarSegregatedOp> m_scalar_segregated_op;

    static_assert(
        ICNS::ndim == AMREX_SPACEDIM,
        "DiffusionOp invoked for scalar PDE type");

    bool use_segregated_op = false;

    DiffusionOp(
        PDEFields& fields, const bool has_overset, const bool mesh_mapping)
    {
        bool use_tensor_op = true;

        amrex::ParmParse pp(fields.field.name() + "_diffusion");
        pp.query("use_tensor_operator", use_tensor_op);
        pp.query("use_segregated_operator", use_segregated_op);

        if (use_tensor_op && use_segregated_op) {
            amrex::Abort(
                "Tensor and segregated operators should not be enabled "
                "simultaneously.");
        }

        if (use_tensor_op) {
            m_tensor_op = std::make_unique<ICNSDiffTensorOp>(
                fields, has_overset, mesh_mapping);
        } else {
            if (use_segregated_op) {
                m_scalar_segregated_op =
                    std::make_unique<ICNSDiffScalarSegregatedOp>(
                        fields, has_overset, mesh_mapping);
            } else {
                m_scalar_op = std::make_unique<ICNSDiffScalarOp>(
                    fields, has_overset, mesh_mapping);
            }
        }
    }

    void compute_diff_term(const FieldState fstate)
    {
        if (m_tensor_op) {
            m_tensor_op->compute_diff_term<Scheme>(fstate);
        } else {
            if (use_segregated_op) {
                m_scalar_segregated_op->compute_diff_term<Scheme>(fstate);
            } else {
                m_scalar_op->compute_diff_term<Scheme>(fstate);
            }
        }
    }

    void linsys_solve(const amrex::Real dt)
    {
        if (m_tensor_op) {
            m_tensor_op->linsys_solve(dt);
        } else {
            if (use_segregated_op) {
                m_scalar_segregated_op->linsys_solve(dt);
            } else {
                m_scalar_op->linsys_solve(dt);
            }
        }
    }
};

} // namespace amr_wind::pde

#endif /* ICNS_DIFFUSION_H */
