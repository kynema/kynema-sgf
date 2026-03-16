#include "RefineCheckpt.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind {
namespace tools {

RefineCheckpt::RefineCheckpt() : incflo() {}

void RefineCheckpt::ErrorEst(
    int lev, amrex::TagBoxArray& tags, amrex::Real, int)
{
    // Tell AMReX to globally refine the mesh
    if (lev < m_orig_ba.size()) {
        tags.setVal(m_orig_ba[lev], amrex::TagBox::SET);
    } else {
        tags.setVal(amrex::TagBox::SET);
    }
}

void RefineCheckpt::read_chkpt_file()
{
    BL_PROFILE("refine-chkpt::read_chkpt_file");
    // Initialize I/O manager to enable restart
    auto& io_mgr = sim().io_manager();
    io_mgr.initialize_io();

    // Ensure that the user has provided a valid checkpoint file
    AMREX_ALWAYS_ASSERT(io_mgr.is_restart());

    ReadCheckpointFile();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Input grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finestLevel());
    }

    m_orig_ba.resize(finestLevel() + 1);
    for (int lev = 0; lev < finestLevel() + 1; ++lev) {
        m_orig_ba[lev] = boxArray(lev);
    }
}

void RefineCheckpt::refine_chkpt_file()
{
    amrex::Print() << "Uniformly refining mesh" << std::endl;
    amrex::Real rstart =
        static_cast<amrex::Real>(amrex::ParallelDescriptor::second());
    regrid(0, sim().time().current_time());
    amrex::Real rend =
        static_cast<amrex::Real>(amrex::ParallelDescriptor::second());
    amrex::Print() << "Refinement time elapsed: " << (rend - rstart)
                   << std::endl;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Refined grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finestLevel());
    }
}

void RefineCheckpt::run_utility()
{
    read_chkpt_file();

    if ((maxLevel() - finestLevel()) != 1) {
        amrex::Print() << "Checkpoint refinement only supported for one "
                          "refinement level at a time\n"
                       << "  Current finest level = " << finestLevel()
                       << "\n  Requested max level  = " << maxLevel()
                       << "\n  Unable to refine checkpoint file" << std::endl;
    } else {
        refine_chkpt_file();
        const int start_level = 1;
        int end_level = finestLevel();
        {
            amrex::ParmParse pp("amr");
            pp.query("max_refined_level", end_level);
        }
        amrex::Print() << "Writing refined levels: " << start_level << " - "
                       << end_level << std::endl;
        sim().io_manager().write_checkpoint_file(start_level, end_level);
    }
}

} // namespace tools
} // namespace amr_wind
