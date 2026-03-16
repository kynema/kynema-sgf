#include <AMReX.H>
#include <AMReX_Print.H>
#include "CheckpointFileUtil.H"
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ctime>

enum format { csv, hdf5 };

void main_main()
{
    std::clock_t time_req;
    const int nargs = amrex::command_argument_count();
    // Define required variables
    std::string checkpointfile;
    format out_format = csv;
    std::string format_str;
    std::string output_file;
    std::string addl_fname;
    amrex::Vector<std::string> addl_fnames{};

    // CLI handling: Collect parameters and help
    for (int i = 1; i <= nargs; i += 2) {
        const std::string& arg = amrex::get_command_argument(i);
        if (arg == "--format" || arg == "-f") {
            // Select format. See enumeration for available options.
            const std::string& arg_format = amrex::get_command_argument(i + 1);
            if (arg_format == "csv") {
                out_format = csv;
                format_str = "csv";
            } else if (arg_format == "hdf5") {
                out_format = hdf5;
                amrex::Abort("HDF5 format has not been implemented yet.");
            } else {
                amrex::Print() << "Invalid format: " << arg_format << std::endl;
                amrex::Abort();
            }
            continue;
        }
        if (arg == "--chk" || arg == "--checkpointfile" || arg == "--file" ||
            arg == "-c") {
            const std::string& chk_tmp = amrex::get_command_argument(i + 1);
            checkpointfile = chk_tmp;
            continue;
        }
        if (arg == "-o" || arg == "--output") {
            output_file = amrex::get_command_argument(i + 1);
        }
        if (arg == "-af" || arg == "--addfield") {
            addl_fname = amrex::get_command_argument(i + 1);
            addl_fnames.push_back(addl_fname);
        }
    }
    amrex::Print() << "Converting checkpoint file " << checkpointfile
                   << " to format " << format_str << std::endl;
    // Help information

    // Load Checkpoint File
    time_req = clock();
    CheckpointFileDataImpl chkptfile(
        checkpointfile, addl_fnames); // Takes care of invalid files

    // Get Meta Data (See Notes 1.)
    int dim = chkptfile.spaceDim();
    int levels = chkptfile.finestLevel() + 1;
    const amrex::Vector<std::string>& var_names =
        chkptfile.varNamesComponents();
    amrex::Real time = chkptfile.time();
    // Least finest level dimensions
    // For info on boxArray, see note 2.
    const amrex::Long nboxes = chkptfile.boxArray(0).size();
    const amrex::Long ncells = chkptfile.boxArray(0).numPts();
    const amrex::Box prob_domain = chkptfile.probDomain(0);
    const auto ncells_domain = prob_domain.d_numPts();
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx = chkptfile.cellSize(0);
    const auto dlo = amrex::lbound(prob_domain);
    const auto dhi = amrex::ubound(prob_domain);
    amrex::IntVect lengths = prob_domain.length(); // hi + 1 for all dims.
    double phis_x = dx[0] * (dhi.x + 1);
    double phis_y = dx[1] * (dhi.y + 1);
    double phis_z = dx[2] * (dhi.z + 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> problo = chkptfile.probLo();
    amrex::Array<amrex::Real, AMREX_SPACEDIM> probhi = chkptfile.probHi();
    // Output Files
    amrex::Print() << "Time: " << time << std::endl
                   << "dx: " << dx << std::endl
                   << "Lengths: " << lengths << std::endl
                   << "lo: " << dlo << " hi: " << dhi << std::endl
                   << "Physical dimensions: (" << phis_x << ", " << phis_y
                   << " , " << phis_z << ")" << std::endl
                   << "Physical coords: " << problo << " to " << probhi
                   << std::endl
                   << "Dimensions: " << dim << std::endl
                   << "Number of Cells: " << ncells << std::endl
                   << "Number of cells/points in domain: " << ncells_domain
                   << std::endl
                   << "Number of Boxes: " << nboxes << std::endl
                   << "levels: " << levels << std::endl;

    if (out_format == csv) {
        // create Header
        /**
         * Header Data in CSV
         *
         * Will follow
         * https://www.w3.org/TR/tabular-data-model/#embedded-metadata
         *
         */
        // Output csv-formatted data
        std::ofstream output_stream;
        output_stream.open(output_file);

        // Create Headers
        output_stream << "x,y,z";
        int num_vars = 0;
        for (auto const& name : var_names) {
            output_stream << "," << name;
            num_vars++;
        }
        output_stream << "\n";
        const amrex::MultiFab& chkptmf = chkptfile.get(0);

        // Loop through MultiFab for data
        for (amrex::MFIter mfi(chkptmf); mfi.isValid(); ++mfi) {
            // Loop through MultiFab object.
            const amrex::Box& bx = mfi.validbox();

            // Seems to be used for masking out non-current box data

            const auto& data = chkptmf.array(mfi); // there is only one box
            // Parallelize the pulling of data
            const amrex::Dim3 lo = amrex::lbound(bx);
            const amrex::Dim3 hi = amrex::ubound(bx);
            for (int z = lo.z; z <= hi.z; ++z) {
                for (int y = lo.y; y <= hi.y; ++y) {
                    // AMREX_PRAGMA_SIMD
                    for (int x = lo.x; x <= hi.x; ++x) {
                        output_stream << x << "," << y << "," << z;
                        for (int n = 0; n < num_vars; n++) {
                            output_stream << "," << data(x, y, z, n);
                        }
                        output_stream << "\n";
                    }
                }
            }
        }
        output_stream.close();
    }

    time_req = clock() - time_req;
    amrex::Print() << "It took " << (float)time_req / CLOCKS_PER_SEC
                   << " seconds to convert." << std::endl;
}

int main(int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}

/*************************/
/**** DEVELOPER NOTES ****/
/*************************/
/*
1. The attributes and member functions are available within theAMReX
    documentation at:
        -
https://amrex-codes.github.io/amrex/doxygen/AMReX__PlotFileUtil_8H_source.html
2. BoxArray contains the information with respect to the boxes/cells.
    The data structure's documentation:
        -
https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1BoxArray.html
3. Documentation on Box object:
        - https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1Box.html
*/
