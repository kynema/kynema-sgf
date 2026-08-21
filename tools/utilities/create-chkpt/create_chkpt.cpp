#include "CreateCheckpt.H"
#include "src/incflo.H"
#include "src/utilities/console_io.H"


int main(int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    using namespace amrex::mpidatatypes;

    kynema_sgf::io::print_banner(MPI_COMM_WORLD, std::cout);

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, []() {
        amrex::ParmParse pp("amrex");
        // Set the defaults so that we throw an exception instead of attempting
        // to generate backtrace files. However, if the user has explicitly set
        // these options in their input files respect those settings.
        if (!pp.contains("throw_exception")) {
            pp.add("throw_exception", 1);
        }
        if (!pp.contains("signal_handling")) {
            pp.add("signal_handling", 0);
        }
    });

    { /* These braces are necessary to ensure amrex::Finalize() can be called
        without explicitly deleting all the incflo member MultiFabs */

        BL_PROFILE("create-chkpt::main");
        kynema_sgf::tools::CreateCheckpt obj;
        obj.run_utility();
    }

    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
