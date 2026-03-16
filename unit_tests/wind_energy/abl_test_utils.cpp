#include "abl_test_utils.H"
#include "aw_test_utils/MeshTest.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

void ABLMeshTest::populate_parameters()
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

        amrex::Vector<amrex::Real> top_velocity{{20.0_rt, 0.0_rt, 0.0_rt}};
        amrex::Vector<amrex::Real> bottom_velocity{{4.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("top_velocity", top_velocity);
        pp.addarr("bottom_velocity", bottom_velocity);
    }

    // Body force
    {
        amrex::ParmParse pp("BodyForce");
        pp.add("type", std::string("oscillatory"));
        amrex::Vector<amrex::Real> source_mag{{1.0_rt, 2.0_rt, 3.0_rt}};
        pp.addarr("magnitude", source_mag);
        pp.add("angular_frequency", 1.0_rt);
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
    }

    // Geostrophic Forcing
    {
        amrex::ParmParse pp("GeostrophicForcing");
        amrex::Vector<amrex::Real> gwind{{10.0_rt, 6.0_rt, 0.0_rt}};
        pp.addarr("geostrophic_wind", gwind);
    }

    // Hurricane Forcing
    {
        amrex::ParmParse pp("HurricaneForcing");
        amrex::Real gradient_wind{40.0_rt};
        amrex::Real radial_distance{40000.0_rt};
        amrex::Real gradient_wind_radial_decay{-0.008_rt};
        amrex::Real gradient_wind_zero_height{18000.};

        pp.add("gradient_wind", gradient_wind);
        pp.add("eyewall_radial_distance", radial_distance);
        pp.add("gradient_wind_radial_decay", gradient_wind_radial_decay);
        pp.add("gradient_wind_zero_height", gradient_wind_zero_height);
    }

    // Rayleigh damping
    {

        amrex::ParmParse pp("RayleighDamping");
        amrex::Real time_scale{40.0_rt};
        amrex::Real length_sloped_damping{200.0_rt};
        amrex::Real length_complete_damping{50.0_rt};
        amrex::Vector<amrex::Real> reference_velocity{
            {12.0_rt, 1.0_rt, -3.0_rt}};
        amrex::Vector<int> fcoord{{1, 0, 1}};

        pp.add("time_scale", time_scale);
        pp.add("length_sloped_damping", length_sloped_damping);
        pp.add("length_complete_damping", length_complete_damping);
        pp.addarr("reference_velocity", reference_velocity);
        pp.addarr("force_coord_directions", fcoord);
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
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);

        pp.add("density", 1.0_rt); // Density
        amrex::Vector<amrex::Real> vel{{20.0_rt, 10.0_rt, 0.0_rt}};
        pp.addarr("velocity", vel);
        amrex::Vector<amrex::Real> grav{{0.0_rt, 0.0_rt, -9.81_rt}};
        pp.addarr("gravity", grav);
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 64}};
        pp.addarr("n_cell", ncell);
    }

    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi{{120.0_rt, 120.0_rt, 1000.0_rt}};
        pp.addarr("prob_hi", probhi);
    }
}

} // namespace amr_wind_tests
