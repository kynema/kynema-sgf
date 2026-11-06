#include "src/core/MLMGOptions.H"

#include "hydro_MacProjector.H"
#include "hydro_NodalProjector.H"

namespace kynema_sgf {

MLMGOptions::MLMGOptions(const std::string& prefix) { parse_options(prefix); }

MLMGOptions::MLMGOptions(
    const std::string& default_prefix, const std::string& custom_prefix)
{
    parse_options(default_prefix);
    parse_options(custom_prefix);
}

void MLMGOptions::parse_options(const std::string& prefix)
{
    amrex::ParmParse pp(prefix);

    // LPInfo options
    pp.query("do_agglomeration", m_lpinfo.do_agglomeration);
    pp.query("do_consolidation", m_lpinfo.do_consolidation);
    pp.query("do_semicoarsening", m_lpinfo.do_semicoarsening);
    pp.query("agg_grid_size", m_lpinfo.agg_grid_size);
    pp.query("con_grid_size", m_lpinfo.con_grid_size);
    pp.query("max_coarsening_level", m_lpinfo.max_coarsening_level);
    pp.query("max_semicoarsening_level", m_lpinfo.max_semicoarsening_level);
    pp.query("deterministic", m_lpinfo.deterministic);

    // LinOp options
    pp.query("max_order", max_order);

    // MLMG options
    m_prefix = prefix;
    std::string multigrid_type;
    if (pp.query("multigrid_type", multigrid_type) != 0) {
        if (multigrid_type == "geometric") {
            m_multigrid_type = amrex::MLMG::MultigridType::geometric;
        } else if (multigrid_type == "algebraic") {
            m_multigrid_type = amrex::MLMG::MultigridType::algebraic;
        } else if (multigrid_type == "hybrid") {
            m_multigrid_type = amrex::MLMG::MultigridType::hybrid;
        } else {
            amrex::Abort(
                "MLMGOptions: invalid " + prefix + ".multigrid_type = " +
                multigrid_type + ". Valid options: geometric algebraic hybrid");
        }
        m_multigrid_type_prefix = prefix;
    }
    pp.query("verbose", m_verbose);
    pp.query("bottom_verbose", m_bottom_verbose);

    pp.query("maxiter", m_max_iter);
    pp.query("mg_rtol", rel_tol);
    pp.query("mg_atol", abs_tol);
    pp.query("fmg_maxiter", m_max_fmg_iters);
    pp.query("num_pre_smooth", m_num_pre_smooth);
    pp.query("num_post_smooth", m_num_post_smooth);
    pp.query("num_final_smooth", m_num_final_smooth);
    pp.query("num_bottom_smooth", m_num_bottom_smooth);

    pp.query("do_fixed_iters", m_do_fixed_iters);

    pp.query("bottom_maxiter", m_bottom_max_iter);
    pp.query("bottom_rtol", m_bottom_rel_tol);
    if (pp.query("bottom_atol", m_bottom_abs_tol) != 0) {
        m_bottom_abs_tol_set = true;
    }
    pp.query("bottom_solver", m_bottom_solver_type);
    pp.query("hypre_namespace", m_hypre_namespace);
    pp.query("hypre_interface", m_hypre_interface);
    pp.query("do_nsolve", m_do_nsolve);
    pp.query("nsolve_grid_size", m_nsolve_grid_size);
}

void MLMGOptions::operator()(amrex::MLMG& mlmg)
{
    mlmg.setVerbose(m_verbose);
    mlmg.setMaxIter(m_max_iter);
    mlmg.setMaxFmgIter(m_max_fmg_iters);

    if (m_do_fixed_iters) {
        mlmg.setFixedIter(m_max_iter);
    }

    mlmg.setNSolve(static_cast<int>(m_do_nsolve));
    mlmg.setNSolveGridSize(m_nsolve_grid_size);
    mlmg.setPreSmooth(m_num_pre_smooth);
    mlmg.setPostSmooth(m_num_post_smooth);
    mlmg.setFinalSmooth(m_num_final_smooth);
    mlmg.setBottomSmooth(m_num_bottom_smooth);

    if (m_multigrid_type != amrex::MLMG::MultigridType::geometric) {
        const std::string key = m_multigrid_type_prefix + ".multigrid_type";
        const auto& linop = mlmg.getLinOp();
        if (!linop.supportsAlgMG() || (linop.getNComp() != 1)) {
            amrex::Abort(
                "MLMGOptions: " + key +
                " = algebraic or hybrid requires a single-component operator "
                "that supports AlgMG, which the " +
                m_prefix + " solve does not have. Set " + m_prefix +
                ".multigrid_type = geometric");
        }
        if (m_do_nsolve) {
            amrex::Abort(
                "MLMGOptions: " + key + " = algebraic or hybrid requires " +
                m_prefix + ".do_nsolve = false");
        }
        if ((m_multigrid_type == amrex::MLMG::MultigridType::hybrid) &&
            m_do_fixed_iters) {
            amrex::Abort(
                "MLMGOptions: " + key + " = hybrid requires " + m_prefix +
                ".do_fixed_iters = false");
        }
    }
    mlmg.setMultigridType(m_multigrid_type);

    // AlgMG solves AMR level 0 with the bottom tolerances. The default
    // absolute bottom tolerance can stop it before the MLMG target is met, so
    // the MLMG iterations stall. Disable it unless the user set it.
    const amrex::Real bottom_abs_tol =
        ((m_multigrid_type != amrex::MLMG::MultigridType::geometric) &&
         !m_bottom_abs_tol_set)
            ? -1.0_rt
            : m_bottom_abs_tol;

    mlmg.setBottomVerbose(m_bottom_verbose);
    mlmg.setBottomTolerance(m_bottom_rel_tol);
    mlmg.setBottomToleranceAbs(bottom_abs_tol);

    if (m_bottom_solver_type == "smoother") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
    } else if (m_bottom_solver_type == "bicg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::bicgstab);
    } else if (m_bottom_solver_type == "cg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::cg);
    } else if (m_bottom_solver_type == "bicgcg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::bicgcg);
    } else if (m_bottom_solver_type == "cgbicg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::cgbicg);
    } else if (m_bottom_solver_type == "amg") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::algmg);
    } else if (m_bottom_solver_type == "hypre") {
#ifdef AMREX_USE_HYPRE
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::hypre);

        mlmg.setHypreOptionsNamespace(m_hypre_namespace);
        if (m_hypre_interface == "ij") {
            mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
        } else if (m_hypre_interface == "semi_structured") {
            mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
        } else if (m_hypre_interface == "structured") {
            mlmg.setHypreInterface(amrex::Hypre::Interface::structed);
        } else {
            amrex::Abort(
                "Invalid hypre interface. Valid options: ij semi_structured "
                "structured");
        }
#else
        amrex::Abort("Kynema-SGF was not built with hypre support");
#endif
    }
}

void MLMGOptions::operator()(Hydro::MacProjector& mac_proj)
{
    operator()(mac_proj.getMLMG());
}

void MLMGOptions::operator()(Hydro::NodalProjector& nodal_proj)
{
    // Only IJ interface supported for NodalProjector
    m_hypre_interface = "ij";

    if ((m_lpinfo.max_coarsening_level > -1) &&
        (m_bottom_solver_type == "hypre")) {
        nodal_proj.getLinOp().setCoarseningStrategy(
            amrex::MLNodeLinOp::CoarseningStrategy::RAP);
    }

    nodal_proj.setVerbose(m_verbose);
    operator()(nodal_proj.getMLMG());
}

} // namespace kynema_sgf
