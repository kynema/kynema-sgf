#include "CoarsenCheckpt.H"
#include "src/utilities/console_io.H"

int main(int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    using namespace amrex::mpidatatypes;

    amr_wind::io::print_banner(MPI_COMM_WORLD, std::cout);

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, []() {
        {
            amrex::ParmParse pp("amrex");
            // Set the defaults so that we throw an exception instead of
            // attempting to generate backtrace files. However, if the user has
            // explicitly set these options in their input files respect those
            // settings.
            if (!pp.contains("throw_exception")) pp.add("throw_exception", 1);
            if (!pp.contains("signal_handling")) pp.add("signal_handling", 0);
        }
        {
            // Preprocess input file
            amrex::ParmParse pp("amr");
            int max_lev_input = 0;
            pp.get("max_level", max_lev_input);
            pp.add("max_level", max_lev_input + 1);
            // Getting the right number of levels is important for geometry init

            amrex::Vector<int> tncell{{0, 0, 0}};
            amrex::Vector<int> ncell{{0, 0, 0}};
            pp.getarr("target_n_cell", tncell);
            pp.getarr("n_cell", ncell);
            // Check n_cell for problems, assume blocking factor of 8
            if ((ncell[0] % 2) % 8 != 0 || (ncell[1] % 2) % 8 != 0 ||
                (ncell[2] % 2) % 8 != 0) {
                amrex::Abort(
                    "Too few cells to coarsen with default blocking factor");
            }
            // Make sure that target matches
            if ((ncell[0] / 2) != tncell[0] || (ncell[1] / 2) != tncell[1] ||
                (ncell[2] / 2) != tncell[2]) {
                amrex::Abort("Coarsened grid does not match target grid");
            }

            // Divide by 2
            ncell[0] /= 2;
            ncell[1] /= 2;
            ncell[2] /= 2;
            pp.addarr("n_cell", ncell);
        }
    });

    {
        BL_PROFILE("coarsen-chkpt::main");
        amr_wind::tools::CoarsenCheckpt obj;
        obj.run_utility();
    }

    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
