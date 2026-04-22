#include <chrono>
#include <ctime>
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/AMRWindVersion.H"
#include "AMReX.H"
#include "AMReX_OpenMP.H"
#include "amr-wind/CFDSim.H"

#ifdef AMR_WIND_USE_NETCDF
#include "netcdf.h"
#ifdef NC_HAVE_META_H
#include "netcdf_meta.h"
#endif
#endif
#ifdef AMR_WIND_USE_MASA
#include "masa.h"
#endif
#ifdef AMREX_USE_HYPRE
#include "HYPRE_config.h"
#endif
#ifdef AMR_WIND_USE_ASCENT
#include "ascent_config.h"
#include "AMReX_REAL.H"

using namespace amrex::literals;
#endif
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace amrex {
const char* buildInfoGetBuildDate();
const char* buildInfoGetComp();
const char* buildInfoGetGitHash(int i);
const char* buildInfoGetCompVersion();
} // namespace amrex

namespace amr_wind::io {

namespace {
const std::string dbl_line = std::string(78, '=') + "\n";
const std::string dash_line = "\n" + std::string(78, '-') + "\n";
} // namespace

void print_usage(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    out << R"doc(Usage:
    amr_wind <input_file> [param=value] [param=value] ...

Required:
    input_file   : Input file with simulation settings

Optional:
    param=value  : Overrides for parameters during runtime
)doc" << '\n';
}

void print_error(MPI_Comm comm, const std::string& msg)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    std::cout << "ERROR: " << msg << '\n';
}

void print_banner(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    auto etime = std::chrono::system_clock::now();
    auto etimet = std::chrono::system_clock::to_time_t(etime);
    amrex::Array<char, 64> time_buf;
    ctime_r(&etimet, time_buf.begin());
    const std::string tstamp(time_buf.begin());
    const std::string dirty_tag = "-" + version::amr_wind_dirty_repo;
    const std::string awind_version = version::amr_wind_version + dirty_tag;
    const std::string awind_git_sha = version::amr_wind_git_sha + dirty_tag;

    // clang-format off
    out << dbl_line
        << "                AMR-Wind (https://github.com/exawind/amr-wind)"
        << '\n' << '\n'
        << "  AMR-Wind version :: " << awind_version << '\n'
        << "  AMR-Wind Git SHA :: " << awind_git_sha << '\n'
        << "  AMReX version    :: " << amrex::Version() << '\n' << '\n'
        << "  Exec. time       :: " << tstamp
        << "  Build time       :: " << amrex::buildInfoGetBuildDate() << '\n'
        << "  C++ compiler     :: " << amrex::buildInfoGetComp()
        << " " << amrex::buildInfoGetCompVersion() << '\n' << '\n'
        << "  MPI              :: "
#ifdef AMREX_USE_MPI
        << "ON    (Num. ranks = " << num_ranks << ")" << std::endl
#else
        << "OFF " << '\n'
#endif
        << "  GPU              :: "
#ifdef AMREX_USE_GPU
        << "ON    "
#if defined(AMREX_USE_CUDA)
        << "(Backend: CUDA)"
#elif defined(AMREX_USE_HIP)
        << "(Backend: HIP)"
#elif defined(AMREX_USE_SYCL)
        << "(Backend: SYCL)"
#endif
        << std::endl
#else
        << "OFF" << '\n'
#endif
        << "  OpenMP           :: "
#ifdef AMREX_USE_OMP
        << "ON    (max threads = " << amrex::OpenMP::get_max_threads()
        << ")" << std::endl
#else
        << "OFF" << '\n'
#endif
        << '\n';

    print_tpls(out);

    out << "           This software is released under the BSD 3-clause license.           "
        << '\n'
        << " See https://github.com/Exawind/amr-wind/blob/development/LICENSE for details. "
        << dash_line << '\n';
    // clang-format on
}

void print_mlmg_header(const std::string& key)
{
    const int name_width = 26;
    amrex::Print() << "\n" << key << '\n';
    amrex::Print() << "  " << std::setw(name_width) << std::left << "System"
                   << std::setw(6) << std::right << "Iters" << std::setw(22)
                   << std::right << "Initial residual" << std::setw(22)
                   << std::right << "Final residual" << '\n'
                   << "  "
                      "--------------------------------------------------------"
                      "--------------------"
                   << '\n';
}

void print_mlmg_info(const std::string& solve_name, const amrex::MLMG& mlmg)
{
    const int name_width = 26;
    amrex::Print() << "  " << std::setw(name_width) << std::left << solve_name
                   << std::setw(6) << std::right << mlmg.getNumIters()
                   << std::setw(22) << std::right << mlmg.getInitResidual()
                   << std::setw(22) << std::right << mlmg.getFinalResidual()
                   << '\n';
}

void print_tpls(std::ostream& out)
{
    amrex::Vector<std::string> tpls;

#ifdef AMR_WIND_USE_NETCDF
    tpls.push_back(std::string("NetCDF    ") + NC_VERSION);
#endif
#ifdef AMREX_USE_HYPRE
    tpls.push_back(std::string("HYPRE     ") + HYPRE_RELEASE_VERSION);
#endif
#ifdef AMR_WIND_USE_OPENFAST
    std::string of_version =
        TOSTRING(OPENFAST_VERSION_MAJOR) "." TOSTRING(OPENFAST_VERSION_MINOR);
    tpls.push_back(std::string("OpenFAST  ") + of_version);
#endif
#ifdef AMR_WIND_USE_MASA
    tpls.push_back(std::string("MASA      ") + MASA_LIB_VERSION);
#endif
#ifdef AMR_WIND_USE_ASCENT
    tpls.push_back(std::string("ASCENT    ") + ASCENT_VERSION);
#endif

    if (!tpls.empty()) {
        out << "  Enabled third-party libraries: ";
        for (const auto& val : tpls) {
            out << "\n    " << val;
        }
        out << '\n' << '\n';
    } else {
        out << "  No additional third-party libraries enabled" << '\n' << '\n';
    }
}

void print_nonlinear_residual(
    const CFDSim& sim, ScratchField& vel_diff, const ScratchField& vel_star)
{
    const int nlevels = sim.repo().num_active_levels();
    const auto& mesh = sim.mesh();

    const auto& velocity_new = sim.pde_manager().icns().fields().field;
    const auto& oset_mask = sim.repo().get_int_field("mask_cell");

    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& velnew_arr = velocity_new(lev).const_arrays();
        const auto& velstar_arr = vel_star(lev).const_arrays();
        const auto& veldiff_arr = vel_diff(lev).arrays();
        const auto& levelmask_arr = level_mask.const_arrays();
        const auto& osetmask_arr = oset_mask(lev).const_arrays();

        amrex::ParallelFor(
            velocity_new(lev), amrex::IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                if (osetmask_arr[nbx](i, j, k) == 0) {
                    veldiff_arr[nbx](i, j, k, n) = 0.;
                } else {
                    veldiff_arr[nbx](i, j, k, n) =
                        (velnew_arr[nbx](i, j, k, n) -
                         velstar_arr[nbx](i, j, k, n)) *
                        levelmask_arr[nbx](i, j, k);
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    amrex::Array<amrex::Real, AMREX_SPACEDIM> rms_vel = {0.0_rt};

    for (int lev = 0; lev < nlevels; ++lev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rms_vel[idim] += vel_diff(lev).norm2(idim);
        }
    }

    amrex::Print() << "Norm of change over fixed point iterations in u:  "
                   << rms_vel[0] << ", v: " << rms_vel[1]
                   << ", w: " << rms_vel[2] << '\n';
}

} // namespace amr_wind::io
