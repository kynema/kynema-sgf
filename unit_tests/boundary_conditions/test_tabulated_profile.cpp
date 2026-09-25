#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/core/FieldRepo.H"
#include "src/physics/udfs/TabulatedProfile.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <fstream>

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

//! Write a profile file that a test can point the boundary condition at
void write_profile(const std::string& fname, const std::string& contents)
{
    std::ofstream outfile(fname);
    outfile << contents;
    outfile.close();
}

/** Fill a field using the tabulated profile and return the largest departure
 *  from the expected linear variation with height
 *
 *  The profile operator only needs the cell index to determine the height, so
 *  it can be exercised over the interior of the domain rather than in the
 *  ghost cells alone.
 */
amrex::Real max_error(
    kynema_sgf::Field& field,
    const amrex::Geometry& geom,
    const kynema_sgf::udf::TabulatedProfile& profile,
    const amrex::Orientation ori,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& slope,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& intercept,
    const amrex::Real zground = 0.0_rt)
{
    const int lev = 0;
    const int ncomp = field.num_comp();
    const auto op = profile.device_instance();
    const auto geomdata = geom.data();
    auto& mfab = field(lev);

    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.validbox();
        const auto& arr = mfab.array(mfi);
        amrex::ParallelFor(
            bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                op(amrex::IntVect{i, j, k}, arr, geomdata, 0.0_rt, ori, n, 0,
                   0);
            });
    }

    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    auto error = amrex::ReduceMax(
        mfab, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& arr) -> amrex::Real {
            amrex::Real err = 0.0_rt;
            amrex::Loop(bx, [=, &err](int i, int j, int k) {
                const auto zco = problo[2] + ((k + 0.5_rt) * dx[2]);
                for (int n = 0; n < ncomp; ++n) {
                    // Below the ground the lowest tabulated value is held
                    const auto zex = amrex::max(zco, zground);
                    const auto expected = intercept[n] + (slope[n] * zex);
                    err = amrex::max(err, std::abs(arr(i, j, k, n) - expected));
                }
            });
            return err;
        });
    amrex::ParallelDescriptor::ReduceRealMax(error);
    return error;
}

//! Write a flat grid file whose ground rises linearly across the domain
void write_terrain(
    const std::string& fname,
    const amrex::Real z_at_ylo,
    const amrex::Real z_at_yhi)
{
    std::ofstream outfile(fname);
    const amrex::Vector<amrex::Real> xs{{0.0_rt, 4.0_rt, 8.0_rt}};
    const amrex::Vector<amrex::Real> ys{{0.0_rt, 4.0_rt, 8.0_rt}};
    outfile << "3\n"
               "3\n";
    for (const auto& x : xs) {
        outfile << x << "\n";
    }
    for (const auto& y : ys) {
        outfile << y << "\n";
    }
    // Indexed [i * ny + j], so x varies slowest
    for (int i = 0; i < 3; ++i) {
        for (const auto& y : ys) {
            outfile << z_at_ylo + ((z_at_yhi - z_at_ylo) * y / 8.0_rt) << "\n";
        }
    }
    outfile.close();
}

} // namespace

class TabulatedProfileTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    //! Declare a field and mark the given faces as inflow
    kynema_sgf::Field& inflow_field(
        const std::string& name,
        const int ncomp,
        const amrex::Vector<amrex::Orientation>& inflow_faces)
    {
        auto& frepo = mesh().field_repo();
        auto& fld = frepo.declare_field(name, ncomp, 1, 1);
        fld.setVal(0.0_rt);
        for (const auto& ori : inflow_faces) {
            fld.bc_type()[ori] = BC::mass_inflow;
        }
        return fld;
    }

    // The default mesh is 8 cells over [0, 8], so heights are 0.5 ... 7.5
    const amrex::Real m_tol = 1.0e-12_rt;
    const amrex::Orientation m_xlo{0, amrex::Orientation::low};
    const amrex::Orientation m_ylo{1, amrex::Orientation::low};
};

TEST_F(TabulatedProfileTest, velocity_from_header)
{
    populate_parameters();
    write_profile(
        "tp_header.txt",
        "# z u v T tke\n"
        "0.0  0.0  3.0  300.0  0.0\n"
        "8.0 16.0 -5.0  308.0  0.8\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_header.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    // u = 2z, v = 3 - z, and w has no column so it is zero
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 0.0_rt},
        {0.0_rt, 3.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, temperature_from_headerless_file)
{
    populate_parameters();
    write_profile(
        "tp_plain.txt",
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_plain.txt"));
    initialize_mesh();

    auto& temp = inflow_field("temperature", 1, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(temp);

    // The fourth column of a headerless file is temperature
    const auto err = max_error(
        temp, mesh().Geom(0), profile, m_xlo, {1.0_rt, 0.0_rt, 0.0_rt},
        {300.0_rt, 0.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, tke_column_is_found_by_field_name)
{
    populate_parameters();
    write_profile(
        "tp_tke.txt",
        "0.0  0.0  3.0  300.0  0.0\n"
        "8.0 16.0 -5.0  308.0  0.8\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_tke.txt"));
    initialize_mesh();

    auto& tke = inflow_field("tke", 1, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(tke);

    const auto err = max_error(
        tke, mesh().Geom(0), profile, m_xlo, {0.1_rt, 0.0_rt, 0.0_rt},
        {0.0_rt, 0.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, each_face_keeps_its_own_profile)
{
    populate_parameters();
    write_profile(
        "tp_x.txt",
        "# z u v T\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    write_profile(
        "tp_y.txt",
        "# z u v T\n"
        "0.0  1.0  0.0  300.0\n"
        "8.0  1.0  8.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_x.txt"));
    }
    {
        amrex::ParmParse pp("ylo");
        pp.add("tabulated_profile_file", std::string("tp_y.txt"));
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo, m_ylo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    const auto err_x = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 0.0_rt},
        {0.0_rt, 3.0_rt, 0.0_rt});
    EXPECT_NEAR(err_x, 0.0_rt, m_tol);

    // The same operator returns the other profile on the other face
    const auto err_y = max_error(
        vel, mesh().Geom(0), profile, m_ylo, {0.0_rt, 1.0_rt, 0.0_rt},
        {1.0_rt, 0.0_rt, 0.0_rt});
    EXPECT_NEAR(err_y, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, face_without_a_profile_uses_the_constant)
{
    populate_parameters();
    write_profile(
        "tp_x_only.txt",
        "# z u v T\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    {
        amrex::ParmParse pp("xlo");
        pp.add("tabulated_profile_file", std::string("tp_x_only.txt"));
    }
    {
        amrex::ParmParse pp("ylo");
        amrex::Vector<amrex::Real> uvw{{7.0_rt, 8.0_rt, 9.0_rt}};
        pp.addarr("velocity", uvw);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo, m_ylo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_ylo, {0.0_rt, 0.0_rt, 0.0_rt},
        {7.0_rt, 8.0_rt, 9.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, ambiguous_column_count_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_bad.txt",
        "0.0  0.0  3.0\n"
        "8.0 16.0 -5.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_bad.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, klaxell_without_a_tke_column_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_no_tke.txt",
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_no_tke.txt"));
    }
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", std::string("KLAxell"));
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, non_monotonic_heights_are_rejected)
{
    populate_parameters();
    write_profile(
        "tp_unsorted.txt",
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n"
        "4.0  8.0 -1.0  304.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_unsorted.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, reversing_normal_velocity_needs_inflow_outflow)
{
    populate_parameters();
    // u enters through xlo low down and leaves higher up, as it does under veer
    write_profile(
        "tp_veer.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0  -4.0  0.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_veer.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, reversing_normal_velocity_is_allowed_on_mixed_face)
{
    populate_parameters();
    write_profile(
        "tp_veer_mio.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0  -4.0  0.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_veer_mio.txt"));
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    vel.setVal(0.0_rt);
    vel.bc_type()[m_xlo] = BC::mass_inflow_outflow;

    const kynema_sgf::udf::TabulatedProfile profile(vel);
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {-1.0_rt, 0.0_rt, 0.0_rt},
        {4.0_rt, 0.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, outflow_everywhere_on_an_inflow_face_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_backwards.txt",
        "# z u v T\n"
        "0.0  -4.0  0.0  300.0\n"
        "8.0  -4.0  0.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_backwards.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, zoffset_lifts_the_profile_to_the_ground)
{
    populate_parameters();
    // Tabulated above ground, on a boundary whose ground sits at z = 2
    write_profile(
        "tp_lift.txt",
        "# z u v T\n"
        "0.0   0.0  1.0  300.0\n"
        "8.0  16.0  1.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_lift.txt"));
        pp.add("zoffset", 2.0_rt);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    // u = 2(z - 2), so the shift moves the whole profile up by the ground
    // height; below the ground the lowest tabulated value is held
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, 0.0_rt, 0.0_rt},
        {-4.0_rt, 1.0_rt, 0.0_rt}, 2.0_rt);
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, each_face_can_sit_on_its_own_ground)
{
    populate_parameters();
    write_profile(
        "tp_ground.txt",
        "# z u v T\n"
        "0.0   0.0  1.0  300.0\n"
        "8.0  16.0  1.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_ground.txt"));
        pp.add("zoffset", 2.0_rt);
    }
    {
        amrex::ParmParse pp("ylo");
        pp.add("tabulated_profile_zoffset", 0.0_rt);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo, m_ylo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    const auto err_x = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, 0.0_rt, 0.0_rt},
        {-4.0_rt, 1.0_rt, 0.0_rt}, 2.0_rt);
    EXPECT_NEAR(err_x, 0.0_rt, m_tol);

    // The face that overrides the offset back to zero is unshifted
    const auto err_y = max_error(
        vel, mesh().Geom(0), profile, m_ylo, {2.0_rt, 0.0_rt, 0.0_rt},
        {0.0_rt, 1.0_rt, 0.0_rt});
    EXPECT_NEAR(err_y, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, a_reversal_the_domain_never_reaches_is_allowed)
{
    populate_parameters();
    // u only turns around above z = 40, far above this 8 m tall domain
    write_profile(
        "tp_high_reversal.txt",
        "# z u v T\n"
        "0.0    4.0  0.0  300.0\n"
        "40.0   4.0  0.0  340.0\n"
        "80.0  -4.0  0.0  380.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_high_reversal.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

TEST_F(TabulatedProfileTest, offset_must_match_the_ground_it_stands_on)
{
    populate_parameters();
    write_profile(
        "tp_g1.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    write_terrain("tp_terrain_flat.amrwind", 3.0_rt, 3.0_rt);
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g1.txt"));
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_terrain_flat.amrwind"));
    }
    initialize_mesh();

    // The ground is at 3 but no offset was given
    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, offset_matching_the_ground_is_accepted)
{
    populate_parameters();
    write_profile(
        "tp_g2.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    write_terrain("tp_terrain_flat2.amrwind", 3.0_rt, 3.0_rt);
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g2.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_terrain_flat2.amrwind"));
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

TEST_F(TabulatedProfileTest, ground_varying_along_the_face_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_g3.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    // Rising across the span, so the xlo face does not stand on level ground
    write_terrain("tp_terrain_slope.amrwind", 0.0_rt, 8.0_rt);
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g3.txt"));
        pp.add("zoffset", 4.0_rt);
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_terrain_slope.amrwind"));
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, an_offset_conflicts_with_an_unaligned_interior)
{
    populate_parameters();
    write_profile(
        "tp_g4.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g4.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    {
        // The interior profile is measured from the bottom of the domain
        amrex::ParmParse pp("ABL");
        pp.add("initial_wind_profile", true);
        pp.add("terrain_aligned_profile", false);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_trailing_word_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e1.txt",
        "# z u v T\\n"
        "0.0 1.0 2.0 300.0 junk\\n"
        "8.0 3.0 4.0 308.0 junk\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e1.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_word_where_a_number_belongs_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e2.txt",
        "# z u v T\\n"
        "0.0 abc 2.0 300.0\\n"
        "8.0 3.0 4.0 308.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e2.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_nan_in_the_file_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e3.txt",
        "# z u v T\\n"
        "0.0 nan 2.0 300.0\\n"
        "8.0 3.0 4.0 308.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e3.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, an_infinity_in_the_file_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e4.txt",
        "# z u v T\\n"
        "0.0 inf 2.0 300.0\\n"
        "8.0 3.0 4.0 308.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e4.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, an_unrepresentable_number_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e5.txt",
        "# z u v T\\n"
        "0.0 1e400 2.0 300.0\\n"
        "8.0 3.0 4.0 308.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e5.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, rows_of_different_widths_are_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e6.txt",
        "# z u v T\\n"
        "0.0 1.0 2.0 300.0\\n"
        "8.0 3.0 4.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e6.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_file_with_no_data_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e7.txt",
        "# z u v T\\n"
        "\\n"
        "# nothing follows\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e7.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_single_height_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e8.txt",
        "# z u v T\\n"
        "0.0 1.0 2.0 300.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e8.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_repeated_column_name_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e9.txt",
        "# z u u T\\n"
        "0.0 1.0 2.0 300.0\\n"
        "8.0 3.0 4.0 308.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e9.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

TEST_F(TabulatedProfileTest, a_height_with_no_values_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e10.txt",
        "# z u v T\\n"
        "0.0\\n"
        "8.0\\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e10.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_THROW(kynema_sgf::udf::TabulatedProfile{vel}, amrex::RuntimeError);
}

} // namespace kynema_sgf_tests
