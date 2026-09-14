#include "CreateCheckpt.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/ncutils/nc_interface.H"
#include "netcdf.h"
#include "netcdf_par.h"

namespace kynema_sgf {
namespace tools {

CreateCheckpt::CreateCheckpt() : incflo() {}


void CreateCheckpt::interpolate()
{

}

#ifdef KYNEMA_SGF_USE_NETCDF
void CreateCheckpt::read_fvcom_data()
{
    amrex::Print() << "Reading FVCOM data" << std::endl;

    // Set up a parmParse to look for "Initialization" inputs.
    amrex::ParmParse pp("Initialization");

    // Get the FVCOM file name.
    std::string m_fvcom_filename;
    pp.query("FVCOM_filename", m_fvcom_filename);
    amrex::Print() << "   -reading file: " << m_fvcom_filename << std::endl;

    // Open the FVCOM file using netcdf utils.
    auto ncdata = ncutils::NCFile::open(m_fvcom_filename, NC_NOWRITE);
    amrex::Print() << "   -file opened..." << std::endl;


    size_t n_elems = ncdata.dim("nele").len();
    size_t n_nodes = ncdata.dim("node").len();
    amrex::Print() << "   -number of elements: " << n_elems << std::endl;
    amrex::Print() << "   -number of nodes:    " << n_nodes << std::endl;

    auto xc = ncdata.var("xc");
    auto yc = ncdata.var("yc");
    auto u = ncdata.var("u");
    auto v = ncdata.var("v");

    amrex::Print() << "    -reading x, y, u, v data..." << std::endl;
    amrex::Vector<amrex::Real> tmp(n_elems, 0.0_rt);
    xc.get(tmp.data(), {0}, {n_elems});
    yc.get(tmp.data(), {0}, {n_elems});
    u.get(tmp.data(), {0}, {n_elems});
    v.get(tmp.data(), {0}, {n_elems});
    amrex::Print() << "   -read data..." << std::endl;


    ncdata.close();
}
#endif

void CreateCheckpt::run_utility()
{
    // Default constructor. Note inheritance: incflo : AmrCore : AmrMesh.
    incflo my_incflo;

    // Initialize data, parameters, arrays and derived internals
    my_incflo.InitData();

    // Set up a parmParse to look for "Initialization" inputs.
    amrex::ParmParse pp("Initialization");
    std::string m_initial_data_type;
    pp.query("data_type", m_initial_data_type);

    amrex::Print() << "data_type = " << m_initial_data_type << std::endl;

    // Call the data reader.
    if ((m_initial_data_type == "fvcom") ||
        (m_initial_data_type == "FVCOM")) {
        read_fvcom_data();
    }
}

} // namespace tools
} // namespace kynema_sgf
