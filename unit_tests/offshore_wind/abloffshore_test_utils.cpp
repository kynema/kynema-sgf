#include "abloffshore_test_utils.H"
#include "aw_test_utils/MeshTest.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

void ABLOffshoreMeshTest::populate_parameters()
{
    // Must populate with default parameters first to signify parameters do not
    // need to be reread within initialize_mesh
    MeshTest::populate_parameters();
    // Initial conditions (Temperature)
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<amrex::Real> theights{
            {0.0_rt, 650.0_rt, 750.0_rt, 1000.0_rt}};
        amrex::Vector<amrex::Real> tvalues{
            {300.0_rt, 300.0_rt, 308.0_rt, 308.75_rt}};
        pp.addarr("temperature_heights", theights);
        pp.addarr("temperature_values", tvalues);
        pp.add("perturb_ref_height", 50.0_rt);
        pp.add("kappa", 0.41_rt);
        pp.add("surface_roughness_z0", 0.1_rt);
    }

    // Initial conditions (Linear profile)
    {
        amrex::ParmParse pp("ABL");
        bool linear_profile = false;
        pp.add("linear_profile", static_cast<int>(linear_profile));

        amrex::Vector<amrex::Real> top_velocity{{20., 0.0_rt, 0.0_rt}};
        amrex::Vector<amrex::Real> bottom_velocity{{4.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("top_velocity", top_velocity);
        pp.addarr("bottom_velocity", bottom_velocity);
    }

    // Transport
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", 300.0_rt);
    }

    // ABL Forcing
    {
        amrex::ParmParse pp("ABLForcing");
        pp.add("abl_forcing_height", 90.0_rt);
        pp.add("abl_forcing_off_height", 10.0_rt);
        pp.add("abl_forcing_ramp_height", 30.0_rt);
    }

    // Geostrophic Forcing
    {
        amrex::ParmParse pp("GeostrophicForcing");
        amrex::Vector<amrex::Real> gwind{{10.0_rt, 6.0_rt, 0.0_rt}};
        pp.addarr("geostrophic_wind", gwind);
        pp.add("wind_forcing_off_height", 10.0_rt);
        pp.add("wind_forcing_ramp_height", 30.0_rt);
    }

    // Coriolis term
    {
        amrex::ParmParse pp("CoriolisForcing");
        pp.add("latitude", 45.0_rt);
    }

    // Needed for initial conditions
    {
        amrex::ParmParse pp("incflo");
        pp.add("use_godunov", 1);
        amrex::Vector<std::string> physics{"MultiPhase", "ABL"};
        pp.addarr("physics", physics);

        pp.add("density", 1.0_rt); // Density
        amrex::Vector<amrex::Real> vel{{20.0_rt, 10.0_rt, 0.0_rt}};
        pp.addarr("velocity", vel);
        amrex::Vector<amrex::Real> grav{{0.0_rt, 0.0_rt, -9.81_rt}};
        pp.addarr("gravity", grav);
    }

    // Water level
    {
        amrex::ParmParse pp("MultiPhase");
        // half cell height
        pp.add("water_level", 0.5_rt * 1000.0_rt / 128.0_rt);
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    // with a negative portion as well for the water to occupy
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 128}};
        pp.addarr("n_cell", ncell);
    }

    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi{{120.0_rt, 120.0_rt, 500.0_rt}};
        pp.addarr("prob_hi", probhi);
        amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, -500.0_rt}};
        pp.addarr("prob_lo", problo);
    }
}

} // namespace amr_wind_tests
