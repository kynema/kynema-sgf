#include "abl_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/wind_energy/ABLFieldInitFile.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
namespace {
void write_ncf()
{
    ncutils::NCFile ncf = ncutils::NCFile::create("abl.nc");
    // Set up dimensions to correspond to ABL test
    ncf.def_dim("nx", 8);
    ncf.def_dim("ny", 8);
    ncf.def_dim("nz", 64);
    // Define vectors that ABLInit will read
    const std::vector<std::string> three_dim{"nx", "ny", "nz"};
    auto uvel = ncf.def_var("uvel", NC_DOUBLE, three_dim);
    auto vvel = ncf.def_var("vvel", NC_DOUBLE, three_dim);
    auto wvel = ncf.def_var("wvel", NC_DOUBLE, three_dim);
    // Populate std vectors
    const std::vector<size_t> start{0, 0, 0};
    const std::vector<size_t> count{8, 8, 64};
    const std::vector<amrex::Real> fill_u(
        static_cast<size_t>(8 * 8 * 64), 0.0_rt);
    const std::vector<amrex::Real> fill_v(
        static_cast<size_t>(8 * 8 * 64), 20.0_rt);
    const std::vector<amrex::Real> fill_w(
        static_cast<size_t>(8 * 8 * 64), 10.0_rt);
    // Populate NetCDF vectors
    uvel.put(fill_u.data(), start, count);
    vvel.put(fill_v.data(), start, count);
    wvel.put(fill_w.data(), start, count);
}
} // namespace

TEST_F(ABLMeshTest, abl_init_netcdf)
{
    populate_parameters();

    // Supply name of input NetCDF file to parser
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_condition_input_file", (std::string) "abl.nc");
    }

    // Create NetCDF file for test
    write_ncf();

    // Set up mock simulation init, starting with mesh and fields
    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto velocity = velocityf.vec_ptrs();

    amr_wind::ABLFieldInitFile ablinitfile;
    run_algorithm(
        mesh().num_levels(), velocity,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            ablinitfile(bx, mesh().Geom(lev), vel, lev);
        });

    const int nlevels = mesh().num_levels();
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Test velocity
    {
        amrex::Vector<amrex::Real> min_vel(3), max_vel(3);
        utils::field_minmax(nlevels, velocity, min_vel, max_vel);
        EXPECT_NEAR(min_vel[0], 0.0_rt, tol);
        EXPECT_NEAR(min_vel[1], 20.0_rt, tol);
        EXPECT_NEAR(min_vel[2], 10.0_rt, tol);
        EXPECT_NEAR(max_vel[0], 0.0_rt, tol);
        EXPECT_NEAR(max_vel[1], 20.0_rt, tol);
        EXPECT_NEAR(max_vel[2], 10.0_rt, tol);
    }
}

TEST_F(ABLMeshTest, abl_init_netcdf_multilevel)
{
    populate_parameters();

    // Supply name of input NetCDF file to parser
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_condition_input_file", (std::string) "abl.nc");
    }
    // If not present, write abl.nc
    std::ifstream f((std::string) "abl.nc");
    if (!f.good()) {
        write_ncf();
    }
    // Set up meshing for multiple levels
    {
        amrex::ParmParse pp("amr");
        pp.add("max_level", 2);
        pp.add("blocking_factor", 2);

        // From abl_test_utils.cpp
        // Grid size is 8 x 8 x 64
        // Grid dimensions are (0.,120.) x (0.,120.) x (0.,1000.)
        std::stringstream ss;
        ss << "2 // Number of levels" << std::endl;
        ss << "1 // Number of boxes at this level" << std::endl;
        ss << "0 0 0 120 120 500" << std::endl;
        ss << "1 // Number of boxes at this level" << std::endl;
        ss << "0 0 0 120 120 200" << std::endl;

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
            new amr_wind::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    // Set up mock simulation init, starting with mesh and fields
    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto velocity = velocityf.vec_ptrs();
    // Default operation needed for fillpatch
    velocityf.set_default_fillpatch_bc(sim().time());

    amr_wind::ABLFieldInitFile ablinitfile;
    bool interp_fine_levels = false;
    const int nlevels = mesh().num_levels();

    for (int lev = 0; lev < nlevels; ++lev) {

        // Fill base level using input file
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(velocityf(lev)); mfi.isValid(); ++mfi) {
            auto vel = velocity[lev]->array(mfi);
            const auto& bx = mfi.validbox();
            interp_fine_levels = ablinitfile(bx, mesh().Geom(lev), vel, lev);
        };

        // Fill the finer levels using coarse data
        if (interp_fine_levels) {
            velocityf.fillpatch_from_coarse(lev, 0.0_rt, velocityf(lev), 0);
        }
    }

    // Test velocity
    {
        const amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
        amrex::Vector<amrex::Real> min_vel(3), max_vel(3);
        utils::field_minmax(nlevels, velocity, min_vel, max_vel);
        EXPECT_NEAR(min_vel[0], 0.0_rt, tol);
        EXPECT_NEAR(min_vel[1], 20.0_rt, tol);
        EXPECT_NEAR(min_vel[2], 10.0_rt, tol);
        EXPECT_NEAR(max_vel[0], 0.0_rt, tol);
        EXPECT_NEAR(max_vel[1], 20.0_rt, tol);
        EXPECT_NEAR(max_vel[2], 10.0_rt, tol);
    }
}

// Clean up ABL NetCDF file if needed
TEST_F(ABLMeshTest, abl_netcdf_cleanup)
{
    // If it exists, remove file
    std::ifstream f((std::string) "abl.nc");
    if (f.good()) {
        remove("abl.nc");
    }
    // Check that file is removed
    std::ifstream ff((std::string) "abl.nc");
    EXPECT_FALSE(ff.good());
}

} // namespace amr_wind_tests
