#include "amr-wind/incflo.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

void incflo::ReadCheckpointFile()
{
    BL_PROFILE("amr-wind::incflo::ReadCheckpointFile()");

    const std::string& restart_file = m_sim.io_manager().restart_file();
    amrex::Print() << "Restarting from checkpoint " << restart_file << '\n';

    amrex::RealArray prob_lo = {0.0_rt};
    amrex::RealArray prob_hi = {0.0_rt};

    /***************************************************************************
     * Load header: set up problem domain (including BoxArray)                 *
     *              allocate incflo memory (incflo::AllocateArrays)            *
     *              (by calling MakeNewLevelFromScratch)
     ***************************************************************************/

    std::string File(restart_file + "/Header");

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // Start reading from checkpoint file

    // Title line
    std::getline(is, line);

    // Finest level
    is >> finest_level;
    amr_wind::ioutils::goto_next_line(is);

    int nstep;
    amrex::Real cur_time, dt_restart;
    // Step count
    is >> nstep;
    amr_wind::ioutils::goto_next_line(is);

    // Current time
    is >> cur_time;
    amr_wind::ioutils::goto_next_line(is);

    m_time.set_restart_time(nstep, cur_time);

    // Time step size
    is >> dt_restart;
    amr_wind::ioutils::goto_next_line(is);

    is >> m_time.delta_t_nm1();
    amr_wind::ioutils::goto_next_line(is);

    is >> m_time.delta_t_nm2();
    amr_wind::ioutils::goto_next_line(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            prob_lo[i++] = static_cast<amrex::Real>(std::stod(word));
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            prob_hi[i++] = static_cast<amrex::Real>(std::stod(word));
        }
    }

    amrex::Vector<amrex::Real> prob_lo_input(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> prob_hi_input(AMREX_SPACEDIM);

    {
        amrex::ParmParse pp("geometry");
        pp.getarr("prob_lo", prob_lo_input);
        pp.getarr("prob_hi", prob_hi_input);
    }

    amrex::Vector<int> n_cell_input(AMREX_SPACEDIM);

    {
        amrex::ParmParse pp("amr");
        pp.getarr("n_cell", n_cell_input);
    }

    amrex::IntVect rep(1, 1, 1);
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        AMREX_ALWAYS_ASSERT(prob_hi[d] > prob_lo[d]);

        const amrex::Real domain_ratio =
            (prob_hi_input[d] - prob_lo_input[d]) / (prob_hi[d] - prob_lo[d]);

        rep[d] = static_cast<int>(domain_ratio);

        constexpr amrex::Real domain_eps = 1.0e-6_rt;
        if (std::abs(static_cast<amrex::Real>(rep[d]) - domain_ratio) >
            domain_eps) {
            amrex::Abort(
                "Domain size changed which indicates replication but there is "
                "either some precision issues or a non integer domain size "
                "change was input");
        }
    }

    bool replicate = (rep != amrex::IntVect::TheUnitVector());

    if (replicate) {
        amrex::Print() << "replicating restart file: " << rep << '\n';
    }

    // Set up problem domain
    amrex::RealBox rb(prob_lo_input.data(), prob_hi_input.data());
    amrex::Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(
            lev, amrex::Geometry(
                     Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                     Geom(lev).isPeriodic()));
    }

    amrex::Vector<amrex::BoxArray> ba_inp(finest_level + 1);
    amrex::Vector<amrex::DistributionMapping> dm_inp(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // read in level 'lev' BoxArray from Header
        ba_inp[lev].readFrom(is);
        amr_wind::ioutils::goto_next_line(is);
    }

    // always use level 0 to check domain size
    constexpr int lev0{0};
    amrex::Box orig_domain(ba_inp[lev0].minimalBox());

    if (replicate) {
        amrex::Print() << " OLD BA had " << ba_inp[lev0].size() << " GRIDS "
                       << '\n';
        amrex::Print() << " OLD Domain" << orig_domain << '\n';
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::BoxList bl;
        for (int k = 0; k < rep[2]; k++) {
            for (int j = 0; j < rep[1]; j++) {
                for (int i = 0; i < rep[0]; i++) {
                    for (int nb = 0; nb < ba_inp[lev].size(); nb++) {
                        amrex::Box b(ba_inp[lev][nb]);
                        amrex::IntVect shift_vec(
                            i * orig_domain.length(0),
                            j * orig_domain.length(1),
                            k * orig_domain.length(2));

                        // equivalent to 2^lev
                        shift_vec *= (1 << lev);

                        b.shift(shift_vec);
                        bl.push_back(b);
                    }
                }
            }
        }
        amrex::BoxArray ba_rep;
        ba_rep.define(bl);

        if (replicate && lev == 0) {

            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                auto new_domain = ba_rep.minimalBox();
                const auto* hi_vect = new_domain.hiVect();

                if (hi_vect[d] + 1 != n_cell_input[d]) {
                    amrex::Abort(
                        "input file error, domain size changed which indicated "
                        "replication, but the amr.n_cell is inconsistent with "
                        "that change in domain size, please adjust amr.n_cell, "
                        "or geometry.prob_lo, or geometry.prob_hi");
                }
            }

            amrex::Print() << " NEW BA had " << ba_rep.size() << " GRIDS "
                           << '\n';
            amrex::Print() << " NEW Domain" << ba_rep.minimalBox() << '\n';
        }

        // Create distribution mapping
        dm_inp[lev].define(ba_inp[lev], amrex::ParallelDescriptor::NProcs());

        amrex::BoxArray ba(ba_rep.simplified());
        ba.maxSize(maxGridSize(lev));
        if (refine_grid_layout) {
            ChopGrids(lev, ba, amrex::ParallelDescriptor::NProcs());
        }
        amrex::DistributionMapping dm =
            amrex::DistributionMapping{ba, amrex::ParallelDescriptor::NProcs()};

        MakeNewLevelFromScratch(lev, m_time.current_time(), ba, dm);
    }

    m_sim.io_manager().read_checkpoint_fields(
        restart_file, ba_inp, dm_inp, rep);
}
