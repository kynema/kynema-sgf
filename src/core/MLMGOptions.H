#ifndef MLMGOPTIONS_H
#define MLMGOPTIONS_H

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MLLinOp.H"
#include "AMReX_MLMG.H"

using namespace amrex::literals;

namespace Hydro {
class NodalProjector;
class MacProjector;
} // namespace Hydro

namespace amr_wind {

/** Interface to control the behavior of AMReX LinearSolvers
 *
 *  MLMGOptions provides a unified interface to set options to the linear
 *  operator object, MLMG solver, external solvers (e.g., hypre), as well as
 *  MacProjector and NodalProjector objects.
 */
struct MLMGOptions
{
    /** Parse user options with a given prefix
     *
     *  \param prefix Prefix used to parse user inputs, e.g., mac_proj
     */
    explicit MLMGOptions(const std::string& prefix);

    /** Parse options in two stages
     *
     *  This constructor parses the default options from `default_prefix` and
     *  then overrides certain options based on a custom prefix.
     *
     *  \param default_prefix Namespace for default options
     *  \param custom_prefix Namespace for override options
     */
    MLMGOptions(
        const std::string& default_prefix, const std::string& custom_prefix);

    void operator()(amrex::MLMG& /*mlmg*/);
    void operator()(Hydro::NodalProjector& /*nodal_proj*/);
    void operator()(Hydro::MacProjector& /*mac_proj*/);

    //! Linear operator options during construction
    amrex::LPInfo& lpinfo() { return m_lpinfo; }

    // Linear operator options
    int max_order{2};

#ifdef BL_USE_FLOAT
    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real rel_tol{std::numeric_limits<amrex::Real>::epsilon() * 1.0e2_rt};

    //! Absolute tolerance for convergence checks
    amrex::Real abs_tol{std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt};
#else
    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real rel_tol{std::numeric_limits<amrex::Real>::epsilon() * 1.0e5_rt};

    //! Absolute tolerance for convergence checks
    amrex::Real abs_tol{std::numeric_limits<amrex::Real>::epsilon() * 1.0e2_rt};
#endif

private:
    void parse_options(const std::string& /*prefix*/);

    //! Linear operator info object
    amrex::LPInfo m_lpinfo;

    // MLMG options

    //! Bottom solver type bicgstab, cg, hypre, etc.
    std::string m_bottom_solver_type{"bicgstab"};

    //! Options namespace for hypre solvers
    std::string m_hypre_namespace{"hypre"};

    //! Options of hypre interface
    std::string m_hypre_interface{"ij"};

#ifdef BL_USE_FLOAT
    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real m_bottom_rel_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e3_rt};

    //! Absolute tolerance for convergence checks
    amrex::Real m_bottom_abs_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e2_rt};
#else
    //! Relative tolerance for convergence of MLMG solvers
    amrex::Real m_bottom_rel_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e12_rt};

    //! Absolute tolerance for convergence checks
    amrex::Real m_bottom_abs_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
#endif

    /** Verbosity of the MLMG and Projector interfaces
     */
    int m_verbose{0};
    int m_max_iter{200};
    int m_max_fmg_iters{0};
    int m_num_pre_smooth{2};
    int m_num_post_smooth{2};
    int m_num_final_smooth{8};
    int m_num_bottom_smooth{0};

    /** Set verbosity for AMReX bottom solvers
     *
     *  When using hypre as the bottom solver, this option is not passed to
     *  hypre. Instead user must set `hypre.verbose` to control the hypre solver
     *  verbosity and `hypre.bamg_verbose` to control BoomerAMG verbosity (when
     *  used as a preconditioner).
     */
    int m_bottom_verbose{0};

    //! Maximum iterations for the bottom solver
    int m_bottom_max_iter{200};

    //! If set to true, AMReX won't abort if linear solver doesn't converge
    //! within max iterations
    bool m_do_fixed_iters{true};

    //! option to use the nsolve bottom solver
    //! added to help with overset masking
    bool m_do_nsolve{false};

    //! set grid size for nsolve
#ifdef AMREX_USE_GPU
    int m_nsolve_grid_size{64};
#else
    int m_nsolve_grid_size{16};
#endif
};

} // namespace amr_wind

#endif /* MLMGOPTIONS_H */
