#include "CoarsenCheckpt.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_PlotFileUtil.H"
#include "amr-wind/utilities/IOManager.H"

namespace {
void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace

namespace amr_wind {
namespace tools {

CoarsenCheckpt::CoarsenCheckpt() : incflo() {}

void CoarsenCheckpt::run_utility()
{
    amrex::Print() << "Reading checkpoint file and adding base coarse level"
                   << std::endl;
    coarsen_chkpt_file();
    const int start_level = 0;
    const int end_level = std::min(start_level, finestLevel() - 1);
    amrex::Print() << "Writing coarsened levels: " << start_level << " - "
                   << end_level << std::endl;
    sim().io_manager().write_checkpoint_file(start_level, end_level);
}

void CoarsenCheckpt::coarsen_chkpt_file()
{

    BL_PROFILE("refine-chkpt::read_chkpt_file");
    // Initialize modified io manager
    sim().io_manager().initialize_io();

    // Ensure that the user has provided a valid checkpoint file
    AMREX_ALWAYS_ASSERT(sim().io_manager().is_restart());

    // Read checkpoint file and add level
    read_chkpt_add_baselevel();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Total (new) grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finestLevel());
    }

    // Average down
    amrex::Print() << "Averaging down to fill new coarse level" << std::endl;
    average_down_all_fields();
}

void CoarsenCheckpt::read_chkpt_add_baselevel()
{
    BL_PROFILE("amr-wind::incflo::ReadCheckpointFile()");

    const std::string& restart_file = sim().io_manager().restart_file();
    amrex::Print() << "Coarsening checkpoint " << restart_file << std::endl;

    amrex::Vector<amrex::Real> prob_lo(AMREX_SPACEDIM);
    amrex::Vector<amrex::Real> prob_hi(AMREX_SPACEDIM);

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
    GotoNextLine(is);

    int nstep;
    amrex::Real cur_time, dt_restart;
    // Step count
    is >> nstep;
    GotoNextLine(is);

    // Current time
    is >> cur_time;
    GotoNextLine(is);

    sim().time().set_restart_time(nstep, cur_time);

    // Time step size
    is >> dt_restart;
    GotoNextLine(is);

    is >> sim().time().delta_t_nm1();
    GotoNextLine(is);

    is >> sim().time().delta_t_nm2();
    GotoNextLine(is);

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

    // Redefine number of levels and base cells
    int finest_level_src = finest_level;
    ++finest_level;
    max_level = finest_level;

    // Set up problem domain
    amrex::RealBox rb(prob_lo.data(), prob_hi.data());
    amrex::Geometry::ResetDefaultProbDomain(rb);
    for (int lev = 0; lev <= max_level; ++lev) {
        SetGeometry(
            lev, amrex::Geometry(
                     Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                     Geom(lev).isPeriodic()));
    }

    amrex::Vector<amrex::BoxArray> ba_inp(finest_level_src + 1);
    amrex::Vector<amrex::DistributionMapping> dm_inp(finest_level_src + 1);

    for (int lev = 0; lev <= finest_level_src; ++lev) {
        // read in level 'lev' BoxArray from Header
        ba_inp[lev].readFrom(is);
        GotoNextLine(is);
    }

    // create base level BoxArray
    const amrex::BoxArray& ba = MakeBaseGrids();
    amrex::DistributionMapping dm =
        amrex::DistributionMapping{ba, amrex::ParallelDescriptor::NProcs()};

    MakeNewLevelFromScratch(0, sim().time().current_time(), ba, dm);

    // Loop and create input levels at higher level numbers
    for (int levsrc = 0; levsrc <= finest_level_src; ++levsrc) {
        int levdst = levsrc + 1;
        // Create distribution mapping
        dm_inp[levsrc].define(
            ba_inp[levsrc], amrex::ParallelDescriptor::NProcs());

        MakeNewLevelFromScratch(
            levdst, sim().time().current_time(), ba_inp[levsrc],
            dm_inp[levsrc]);
    }

    // Read checkpoint fields into non-coarsened destination levels
    amrex::IntVect rep(1, 1, 1); // no replication
    read_checkpoint_fields_offset(restart_file, ba_inp, dm_inp, rep);
}

void CoarsenCheckpt::read_checkpoint_fields_offset(
    const std::string& restart_file,
    const amrex::Vector<amrex::BoxArray>& ba_chk,
    const amrex::Vector<amrex::DistributionMapping>& dm_chk,
    const amrex::IntVect& rep)
{
    BL_PROFILE("amr-wind::IOManager::read_checkpoint_fields");

    // Track set of fields that might be missing at this level
    std::set<std::string> missing;
    const std::string level_prefix = "Level_";
    const int nlevels = sim().mesh().finestLevel() + 1;

    // always use the level 0 domain
    amrex::Box orig_domain(ba_chk[0].minimalBox());

    for (int levsrc = 0; levsrc < nlevels - 1; ++levsrc) {
        const int levdst = levsrc + 1;
        for (auto* fld : sim().io_manager().checkpoint_fields()) {
            auto& field = *fld;
            const auto& fab_file = amrex::MultiFabFileFullPrefix(
                levsrc, restart_file, level_prefix, field.name());

            // Fields might be registered for checkpoint but might not be
            // necessary for actually performing the simulation. Check if the
            // field exists before attempting to read the restart field.
            if (!amrex::VisMF::Exist(fab_file)) {
                missing.insert(field.name());
                continue;
            }

            auto& mfab = field(levdst);
            const auto& ba_fab = amrex::convert(ba_chk[levsrc], mfab.ixType());
            if (mfab.boxArray() == ba_fab &&
                mfab.DistributionMap() == dm_chk[levsrc]) {
                amrex::VisMF::Read(
                    field(levdst),
                    amrex::MultiFabFileFullPrefix(
                        levsrc, restart_file, level_prefix, field.name()));
            } else {
                amrex::MultiFab tmp(
                    ba_fab, dm_chk[levdst], mfab.nComp(), mfab.nGrowVect());
                amrex::VisMF::Read(
                    tmp, amrex::MultiFabFileFullPrefix(
                             levsrc, restart_file, level_prefix, field.name()));

                for (int k = 0; k < rep[2]; k++) {
                    for (int j = 0; j < rep[1]; j++) {
                        for (int i = 0; i < rep[0]; i++) {

                            amrex::IntVect shift_vec(
                                i * orig_domain.length(0),
                                j * orig_domain.length(1),
                                k * orig_domain.length(2));

                            // equivalent to 2^lev
                            shift_vec *= (1 << levdst);

                            tmp.shift(shift_vec);
                            mfab.ParallelCopy(tmp);
                            tmp.shift(-shift_vec);
                        }
                    }
                }

                mfab.setBndry(0.0);
            }
        }
    }

    // If fields were missing, print diagnostic message.
    if (!missing.empty()) {
        amrex::Print() << "\nWARNING: The following fields were missing in the "
                          "restart file for one or more levels. Please check "
                          "your restart file and inputs."
                       << std::endl
                       << "Missing checkpoint fields: " << std::endl;
        for (const auto& ff : missing) {
            amrex::Print() << "  - " << ff << std::endl;
        }
        amrex::Print() << std::endl;
        if (true) {
            // Need to have all fields
            amrex::Abort("Missing fields in restart file.");
        }
    }
}

void CoarsenCheckpt::average_down_all_fields()
{
    // Target level is the new, unpopulated base level
    int lev = 0;
    for (auto* fld : sim().io_manager().checkpoint_fields()) {
        auto& field = *fld;
        amrex::average_down(
            field(lev + 1), field(lev), 0, field.num_comp(),
            sim().mesh().refRatio(lev));
    }
}

} // namespace tools
} // namespace amr_wind
