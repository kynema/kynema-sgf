#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/BCInterface.H"
#include "src/boundary_conditions/scalar_bcs.H"
#include "src/core/FieldBCOps.H"
#include "src/boundary_conditions/velocity_bcs.H"
#include "src/equation_systems/PDEHelpers.H"
#include "src/equation_systems/SchemeTraits.H"
#include "src/equation_systems/tke/TKE.H"
#include "src/core/FieldRepo.H"
#include "src/physics/udfs/TabulatedProfile.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <cstdio>
#include <fstream>
#include <limits>

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

/** Check that building the profile aborts, and for the expected reason
 *
 *  Several checks can reject the same broken file, so a bare EXPECT_THROW
 *  would pass on whichever fires first. Matching the message pins the test to
 *  the check it is named after.
 *
 *  \param field Field whose boundary profile is built
 *  \param message Part of the abort message the check under test gives
 */
void expect_abort_with(
    const kynema_sgf::Field& field, const std::string& message)
{
    try {
        const kynema_sgf::udf::TabulatedProfile profile(field);
        ADD_FAILURE() << "expected an abort mentioning: " << message;
    } catch (const amrex::RuntimeError& err) {
        const std::string what = err.what();
        EXPECT_NE(what.find(message), std::string::npos)
            << "aborted for another reason: " << what;
    }
}

/** Set the normal velocity through xlo to enter below a given height and
 *  leave above it, ghost cells included, as a veering inflow does
 *
 *  \param vel Velocity field to fill
 *  \param kmid First cell index, counted up from the bottom, that flows out
 */
void set_inflow_outflow_velocity(kynema_sgf::Field& vel, const int kmid)
{
    auto& mfab = vel(0);
    for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const auto& arr = mfab.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            arr(i, j, k, 0) = (k < kmid) ? 1.0_rt : -1.0_rt;
            arr(i, j, k, 1) = 0.0_rt;
            arr(i, j, k, 2) = 0.0_rt;
        });
    }
}

/** Largest departure of the xlo ghost cells from the profile where the flow
 *  enters and from the adjacent interior value where it leaves
 *
 *  \param field Scalar field whose ghost cells are checked
 *  \param kmid First cell index, counted up from the bottom, that flows out
 *  \param slope Rate at which the tabulated profile increases with height
 *  \param interior Value held in every interior cell
 */
amrex::Real xlo_ghost_error(
    kynema_sgf::Field& field,
    const int kmid,
    const amrex::Real slope,
    const amrex::Real interior)
{
    const auto& domain = field.repo().mesh().Geom(0).Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    auto error = amrex::ReduceMax(
        field(0), 1,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& arr) -> amrex::Real {
            amrex::Real err = 0.0_rt;
            amrex::Loop(bx, [=, &err](int i, int j, int k) {
                if ((i != dlo.x - 1) || (j < dlo.y) || (j > dhi.y) ||
                    (k < dlo.z) || (k > dhi.z)) {
                    return;
                }
                // The grid is 1 m, so the cell center height is k + 0.5
                const auto expected =
                    (k < kmid) ? (slope * (k + 0.5_rt)) : interior;
                err = amrex::max(err, std::abs(arr(i, j, k) - expected));
            });
            return err;
        });
    amrex::ParallelDescriptor::ReduceRealMax(error);
    return error;
}

/** Value the profile operator puts at one index of one face
 *
 *  \param profile Profile under test
 *  \param geom Geometry of the level
 *  \param ori Face being filled
 *  \param iv Index filled
 *  \param comp Component filled
 *  \param face_dir Face direction of a face-centered field, or -1
 *  \return The value
 */
amrex::Real value_at(
    const kynema_sgf::udf::TabulatedProfile& profile,
    const amrex::Geometry& geom,
    const amrex::Orientation ori,
    const amrex::IntVect& iv,
    const int comp,
    const int face_dir = -1)
{
    const auto op = profile.device_instance(face_dir);
    const auto geomdata = geom.data();
    const amrex::Box bx(iv, iv);
    amrex::FArrayBox fab(bx, 1, amrex::The_Pinned_Arena());
    const auto& arr = fab.array();
    // Component comp of the field, written to component 0 of the box
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        op(amrex::IntVect{i, j, k}, arr, geomdata, 0.0_rt, ori, 0, 0, comp);
    });
    amrex::Gpu::streamSynchronize();
    return fab(iv, 0);
}

/** Fill a field using the tabulated profile and return the largest departure
 *  from the expected linear variation with height
 *
 *  The profile operator only needs the cell index to determine the height, so
 *  it can be exercised over the interior of the domain rather than in the
 *  ghost cells alone. With face_dir = 2 the indices are taken as the z nodes
 *  of a face-centered field, as for w on its z faces.
 */
amrex::Real max_error(
    kynema_sgf::Field& field,
    const amrex::Geometry& geom,
    const kynema_sgf::udf::TabulatedProfile& profile,
    const amrex::Orientation ori,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& slope,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& intercept,
    const amrex::Real zground = 0.0_rt,
    const int face_dir = -1)
{
    const int lev = 0;
    const int ncomp = field.num_comp();
    const auto op = profile.device_instance(face_dir);
    // A face-centered field is evaluated on the nodes of its face direction
    const amrex::Real shift = (face_dir == 2) ? 0.0_rt : 0.5_rt;
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
                const auto zco = problo[2] + ((k + shift) * dx[2]);
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

//! Write a flat grid file whose ground varies along y only, knot by knot
void write_terrain_along_y(
    const std::string& fname,
    const amrex::Vector<amrex::Real>& ys,
    const amrex::Vector<amrex::Real>& zs)
{
    std::ofstream outfile(fname);
    const amrex::Vector<amrex::Real> xs{{0.0_rt, 8.0_rt}};
    outfile << xs.size() << "\n" << ys.size() << "\n";
    for (const auto& x : xs) {
        outfile << x << "\n";
    }
    for (const auto& y : ys) {
        outfile << y << "\n";
    }
    // Indexed [i * ny + j], so x varies slowest
    for (int i = 0; i < xs.size(); ++i) {
        for (const auto& z : zs) {
            outfile << z << "\n";
        }
    }
    outfile.close();
}

//! Write a flat grid file whose ground rises linearly along x only
void write_terrain_along_x(
    const std::string& fname,
    const amrex::Real z_at_xlo,
    const amrex::Real z_at_xhi)
{
    std::ofstream outfile(fname);
    const amrex::Vector<amrex::Real> xs{{0.0_rt, 8.0_rt}};
    const amrex::Vector<amrex::Real> ys{{0.0_rt, 8.0_rt}};
    outfile << "2\n"
               "2\n";
    for (const auto& x : xs) {
        outfile << x << "\n";
    }
    for (const auto& y : ys) {
        outfile << y << "\n";
    }
    // Indexed [i * ny + j], so x varies slowest
    for (const auto& x : xs) {
        for (int j = 0; j < ys.size(); ++j) {
            outfile << z_at_xlo + ((z_at_xhi - z_at_xlo) * x / 8.0_rt) << "\n";
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

    //! Declare a field, mark the given faces as inflow and select the
    //! tabulated profile on them (<face>.<field>.inflow_type)
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
            amrex::ParmParse pp(kynema_sgf::bcnames[ori]);
            pp.add(name + ".inflow_type", std::string("TabulatedProfile"));
        }
        return fld;
    }

    // The default mesh is 8 cells over [0, 8], so heights are 0.5 ... 7.5
    // Relative to machine precision so it holds in single precision too; the
    // profiles are linear, so the interpolation itself is exact
    const amrex::Real m_tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    const amrex::Orientation m_xlo{0, amrex::Orientation::low};
    const amrex::Orientation m_ylo{1, amrex::Orientation::low};
};

// A note that starts with z is not the header, before or after the real one.
// The header swaps u and v, so only the real header gives u = 2z, v = 3 - z
TEST_F(TabulatedProfileTest, a_note_before_the_header_is_not_the_header)
{
    populate_parameters();
    write_profile(
        "tp_note1.txt",
        "# z is the height above ground in meters\n"
        "# z v u T\n"
        "0.0  3.0  0.0  300.0\n"
        "8.0 -5.0 16.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_note1.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 0.0_rt},
        {0.0_rt, 3.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

TEST_F(TabulatedProfileTest, a_note_after_the_header_is_not_the_header)
{
    populate_parameters();
    write_profile(
        "tp_note2.txt",
        "# z v u T\n"
        "# z values are in meters\n"
        "0.0  3.0  0.0  300.0\n"
        "8.0 -5.0 16.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_note2.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 0.0_rt},
        {0.0_rt, 3.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

// Notes that start with z do not make a headerless file need a header: with
// no candidate as wide as the data and none naming a known column, the file
// is read as z u v T
TEST_F(TabulatedProfileTest, a_note_alone_leaves_the_file_headerless)
{
    populate_parameters();
    write_profile(
        "tp_note3.txt",
        "# z is the height above ground\n"
        "# z is the z coordinate\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_note3.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 0.0_rt},
        {0.0_rt, 3.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

// A header that names known columns but does not fit the data is a mistake,
// reported with its line, rather than a note to skip
TEST_F(TabulatedProfileTest, a_header_too_wide_for_the_data_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_note4.txt",
        "# z u v T tke\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_note4.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "tp_note4.txt line 1 reads as a header naming 5 columns, but the data "
        "has 4");
}

// Nor with a header: the RANS readers take none, so sharing is refused
TEST_F(TabulatedProfileTest, a_header_does_not_allow_sharing_the_rans_file)
{
    populate_parameters();
    write_profile(
        "tp_rans3.txt",
        "# z u v T tke\n"
        "0.0 1.0 2.0 300.0 0.5\n"
        "8.0 3.0 4.0 308.0 0.4\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("rans_1dprofile_file", std::string("tp_rans3.txt"));
    }
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_rans3.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "Use a separate file for the inflow profile");
}

// A note alone does not hide that the file is the RANS profile
TEST_F(TabulatedProfileTest, a_note_does_not_hide_the_rans_file)
{
    populate_parameters();
    write_profile(
        "tp_rans2.txt",
        "# z is the height above ground\n"
        "0.0 1.0 2.0 0.0 0.5\n"
        "8.0 3.0 4.0 0.0 0.4\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("rans_1dprofile_file", std::string("tp_rans2.txt"));
    }
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_rans2.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "is also used as ABL.rans_1dprofile_file");
}

TEST_F(TabulatedProfileTest, a_missing_profile_file_is_rejected)
{
    populate_parameters();
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_no_such_profile.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "cannot open profile file tp_no_such_profile.txt");
}

// Header names are matched without regard to case, and T, theta, temp and
// temperature all name the temperature column
TEST_F(TabulatedProfileTest, header_names_ignore_case_and_accept_aliases)
{
    populate_parameters();
    write_profile(
        "tp_alias.txt",
        "# Z Theta U V\n"
        "0.0  300.0  0.0  3.0\n"
        "8.0  308.0 16.0 -5.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_alias.txt"));
    initialize_mesh();

    auto& temp = inflow_field("temperature", 1, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(temp);
    // T = 300 + z comes from the second column, so only the header gives it:
    // read as headerless z u v T, the fourth column would be taken instead
    const auto err = max_error(
        temp, mesh().Geom(0), profile, m_xlo, {1.0_rt, 0.0_rt, 0.0_rt},
        {300.0_rt, 0.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol * 300.0_rt);
}

// '## z ...' and '# z, u, ...' are headers too: leading '#' characters and
// commas are not part of the names. T is in the fifth column, where a
// headerless file would have tke, so only the header gives T = 300 + z
TEST_F(TabulatedProfileTest, headers_with_extra_hashes_or_commas_are_read)
{
    populate_parameters();
    write_profile(
        "tp_hash.txt",
        "## z, u, v, w, T\n"
        "0.0  0.0  3.0  0.0  300.0\n"
        "8.0 16.0 -5.0  0.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_hash.txt"));
    initialize_mesh();

    auto& temp = inflow_field("temperature", 1, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(temp);
    const auto err = max_error(
        temp, mesh().Geom(0), profile, m_xlo, {1.0_rt, 0.0_rt, 0.0_rt},
        {300.0_rt, 0.0_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol * 300.0_rt);
}

// A comment that names the columns but does not start with z is reported
// rather than skipped, which would have the columns guessed
TEST_F(TabulatedProfileTest, a_header_not_starting_with_z_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_height.txt",
        "# height u v w T\n"
        "0.0  0.0  3.0  0.0  300.0\n"
        "8.0 16.0 -5.0  0.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_height.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "tp_height.txt line 1 looks like a header for the 5 columns of "
        "the data, but a header must start with 'z'");
}

// One known column is enough: '# height u humidity speed' over four columns
// would otherwise be read as z u v T, with humidity taken as v
TEST_F(TabulatedProfileTest, a_header_naming_one_known_column_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_height1.txt",
        "# height u humidity speed\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_height1.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "tp_height1.txt line 1 looks like a header for the 4 columns of "
        "the data, but a header must start with 'z'");
}

// A note as wide as the data is taken as the header; the error then names
// its line so that the note can be reworded
TEST_F(TabulatedProfileTest, a_note_taken_as_header_is_named)
{
    populate_parameters();
    write_profile(
        "tp_note5.txt",
        "# z is height above ground\n"
        "0.0  0.0  3.0  300.0  0.1\n"
        "8.0 16.0 -5.0  308.0  0.1\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_note5.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "has no u column, needed for the velocity boundary condition on "
        "xlo (columns named by line 1");
}

// The BC setup refuses different UDFs on the mass_inflow and the
// mass_inflow_outflow faces of a field, since one fill operator serves both
TEST_F(TabulatedProfileTest, different_udfs_on_the_two_face_types_are_rejected)
{
    populate_parameters();
    for (const auto& face : {"ylo", "zlo", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("type", std::string("mass_inflow"));
        pp.add("velocity.inflow_type", std::string("CustomVelocity"));
        pp.add("tracer.inflow_type", std::string("CustomScalar"));
        pp.add("tracer", 0.0_rt);
    }
    {
        // Read by CustomVelocity, so that a missing check fails the test
        // rather than stopping on an input error
        amrex::ParmParse pp("CustomVelocity");
        amrex::Vector<amrex::Real> uvw{{1.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
    }
    {
        amrex::ParmParse pp("xhi");
        pp.add("type", std::string("mass_inflow_outflow"));
        pp.add("velocity.inflow_outflow_type", std::string("TabulatedProfile"));
        pp.add("tracer.inflow_outflow_type", std::string("TabulatedProfile"));
        pp.add("tracer", 0.0_rt);
    }
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    kynema_sgf::BCVelocity vbc(vel);
    vbc();
    try {
        kynema_sgf::vel_bc::register_velocity_dirichlet(
            vel, mesh(), time(), vbc.get_dirichlet_udfs());
        ADD_FAILURE() << "expected the velocity UDFs to be refused";
    } catch (const amrex::RuntimeError& err) {
        EXPECT_NE(
            std::string(err.what())
                .find(
                    "differ; one UDF fills every inflow "
                    "face"),
            std::string::npos)
            << err.what();
    }

    auto& tracer = frepo.declare_field("tracer", 1, 1, 1);
    kynema_sgf::BCScalar sbc(tracer);
    sbc(0.0_rt);
    try {
        kynema_sgf::scalar_bc::register_scalar_dirichlet(
            tracer, mesh(), time(), sbc.get_dirichlet_udfs());
        ADD_FAILURE() << "expected the tracer UDFs to be refused";
    } catch (const amrex::RuntimeError& err) {
        EXPECT_NE(
            std::string(err.what())
                .find(
                    "differ; one UDF fills every inflow "
                    "face"),
            std::string::npos)
            << err.what();
    }
}

namespace {
//! Slip walls on the y and z faces, and the given types on xlo and xhi
void set_x_faces(const std::string& xlo_type, const std::string& xhi_type)
{
    for (const auto& face : {"ylo", "zlo", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    amrex::ParmParse("xlo").add("type", xlo_type);
    amrex::ParmParse("xhi").add("type", xhi_type);
    // Read by CustomVelocity, so that a missing check fails the test rather
    // than stopping on an input error
    amrex::Vector<amrex::Real> uvw{{1.0_rt, 0.0_rt, 0.0_rt}};
    amrex::ParmParse("CustomVelocity").addarr("velocity", uvw);
}

//! Register the inflow UDFs of a field and return the error, if any
template <typename BCType, typename Register>
std::string register_error(
    kynema_sgf::Field& fld, const Register& register_fn, const bool vector)
{
    BCType bc(fld);
    if (vector) {
        bc();
    } else {
        bc(0.0_rt);
    }
    try {
        register_fn(bc.get_dirichlet_udfs());
    } catch (const amrex::RuntimeError& err) {
        return err.what();
    }
    return "";
}
} // namespace

// The one fill operator calls a custom UDF on every inflow face, so a face of
// the other kind that gives a constant is refused rather than overwritten
TEST_F(TabulatedProfileTest, a_custom_udf_next_to_a_constant_face_is_rejected)
{
    populate_parameters();
    set_x_faces("mass_inflow", "mass_inflow_outflow");
    {
        amrex::ParmParse pp("xlo");
        pp.add("velocity.inflow_type", std::string("CustomVelocity"));
        pp.add("tracer.inflow_type", std::string("CustomScalar"));
        pp.add("tracer", 0.0_rt);
    }
    {
        amrex::ParmParse pp("xhi");
        amrex::Vector<amrex::Real> uvw{{1.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
        pp.add("tracer", 0.0_rt);
    }
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    const auto verr = register_error<kynema_sgf::BCVelocity>(
        vel,
        [&](const auto& udfs) {
            kynema_sgf::vel_bc::register_velocity_dirichlet(
                vel, mesh(), time(), udfs);
        },
        true);
    EXPECT_NE(
        verr.find(
            "xhi.velocity.inflow_outflow_type is ConstDirichlet, but "
            "CustomVelocity fills every inflow face"),
        std::string::npos)
        << verr;

    auto& tracer = frepo.declare_field("tracer", 1, 1, 1);
    const auto serr = register_error<kynema_sgf::BCScalar>(
        tracer,
        [&](const auto& udfs) {
            kynema_sgf::scalar_bc::register_scalar_dirichlet(
                tracer, mesh(), time(), udfs);
        },
        false);
    EXPECT_NE(
        serr.find(
            "xhi.tracer.inflow_outflow_type is ConstDirichlet, but "
            "CustomScalar fills every inflow face"),
        std::string::npos)
        << serr;
}

// The same holds between two faces of the same kind
TEST_F(TabulatedProfileTest, a_custom_udf_must_be_selected_on_every_inflow_face)
{
    populate_parameters();
    set_x_faces("mass_inflow", "mass_inflow");
    amrex::ParmParse("xlo").add(
        "velocity.inflow_type", std::string("CustomVelocity"));
    {
        amrex::ParmParse pp("xhi");
        amrex::Vector<amrex::Real> uvw{{-1.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
    }
    initialize_mesh();

    auto& vel = mesh().field_repo().declare_field("velocity", 3, 1, 1);
    const auto verr = register_error<kynema_sgf::BCVelocity>(
        vel,
        [&](const auto& udfs) {
            kynema_sgf::vel_bc::register_velocity_dirichlet(
                vel, mesh(), time(), udfs);
        },
        true);
    EXPECT_NE(
        verr.find(
            "xhi.velocity.inflow_type is ConstDirichlet, but "
            "CustomVelocity fills every inflow face"),
        std::string::npos)
        << verr;
}

// A custom UDF on every inflow face passes the check (CustomVelocity itself
// is a template that stops when built), and TabulatedProfile may sit next to
// a constant face, since it falls back to the constant face by face
TEST_F(TabulatedProfileTest, udfs_that_can_fill_their_faces_are_registered)
{
    populate_parameters();
    set_x_faces("mass_inflow", "pressure_outflow");
    amrex::ParmParse("xlo").add(
        "velocity.inflow_type", std::string("CustomVelocity"));
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    EXPECT_EQ(
        register_error<kynema_sgf::BCVelocity>(
            vel,
            [&](const auto& udfs) {
                kynema_sgf::bc_udf::check_inflow_udf_faces(
                    vel, udfs, "Velocity BC");
            },
            true),
        "");

    write_profile(
        "tp_mixed_const.txt",
        "# z u v T tracer\n"
        "0.0  1.0  0.0  300.0  1.0\n"
        "8.0  2.0  0.0  308.0  2.0\n");
    amrex::ParmParse("TabulatedProfile")
        .add("filename", std::string("tp_mixed_const.txt"));
    amrex::ParmParse("xlo").add(
        "tracer.inflow_type", std::string("TabulatedProfile"));
    amrex::ParmParse("xhi").add("type", std::string("mass_inflow_outflow"));
    amrex::ParmParse("xhi").add("tracer", 0.0_rt);
    auto& tracer = frepo.declare_field("tracer", 1, 1, 1);
    EXPECT_EQ(
        register_error<kynema_sgf::BCScalar>(
            tracer,
            [&](const auto& udfs) {
                kynema_sgf::scalar_bc::register_scalar_dirichlet(
                    tracer, mesh(), time(), udfs);
            },
            false),
        "");
}

// A face that names ConstDirichlet is a face that does not select the
// profile: it keeps its constant next to a TabulatedProfile face, including
// through the BC setup, which gathers the UDFs of all faces
TEST_F(TabulatedProfileTest, an_explicit_constant_face_keeps_its_constant)
{
    populate_parameters();
    set_x_faces("mass_inflow", "mass_inflow");
    write_profile(
        "tp_explicit_const.txt",
        "# z u v T\n"
        "0.0  1.0  0.0  300.0\n"
        "8.0  2.0  0.0  308.0\n");
    amrex::ParmParse("TabulatedProfile")
        .add("filename", std::string("tp_explicit_const.txt"));
    amrex::ParmParse("xlo").add(
        "velocity.inflow_type", std::string("TabulatedProfile"));
    {
        amrex::ParmParse pp("xhi");
        pp.add("velocity.inflow_type", std::string("ConstDirichlet"));
        amrex::Vector<amrex::Real> uvw{{-1.0_rt, 0.5_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
    }
    initialize_mesh();

    auto& vel = mesh().field_repo().declare_field("velocity", 3, 1, 1);
    EXPECT_EQ(
        register_error<kynema_sgf::BCVelocity>(
            vel,
            [&](const auto& udfs) {
                kynema_sgf::vel_bc::register_velocity_dirichlet(
                    vel, mesh(), time(), udfs);
            },
            true),
        "");

    const kynema_sgf::udf::TabulatedProfile profile(vel);
    const amrex::Orientation xhi{0, amrex::Orientation::high};
    const auto err = max_error(
        vel, mesh().Geom(0), profile, xhi, {0.0_rt, 0.0_rt, 0.0_rt},
        {-1.0_rt, 0.5_rt, 0.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

// Next to a UDF that fills every inflow face, a face that names
// ConstDirichlet is refused like one that names nothing
TEST_F(TabulatedProfileTest, an_explicit_constant_face_next_to_a_custom_udf)
{
    populate_parameters();
    set_x_faces("mass_inflow", "mass_inflow");
    amrex::ParmParse("xlo").add(
        "velocity.inflow_type", std::string("CustomVelocity"));
    {
        amrex::ParmParse pp("xhi");
        pp.add("velocity.inflow_type", std::string("ConstDirichlet"));
        amrex::Vector<amrex::Real> uvw{{-1.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
    }
    initialize_mesh();

    auto& vel = mesh().field_repo().declare_field("velocity", 3, 1, 1);
    const auto verr = register_error<kynema_sgf::BCVelocity>(
        vel,
        [&](const auto& udfs) {
            kynema_sgf::vel_bc::register_velocity_dirichlet(
                vel, mesh(), time(), udfs);
        },
        true);
    EXPECT_NE(
        verr.find(
            "xhi.velocity.inflow_type is ConstDirichlet, but "
            "CustomVelocity fills every inflow face"),
        std::string::npos)
        << verr;
}

// The constant density of ICNS registers no inflow UDF, so selecting the
// profile for it is refused in its BC setup rather than ignored
TEST_F(TabulatedProfileTest, a_tabulated_constant_density_is_rejected)
{
    populate_parameters();
    for (const auto& face : {"ylo", "zlo", "xhi", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("type", std::string("mass_inflow"));
        amrex::Vector<amrex::Real> uvw{{1.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
        pp.add("density", 1.0_rt);
        pp.add("density.inflow_type", std::string("TabulatedProfile"));
    }
    initialize_mesh();

    try {
        sim().pde_manager().register_icns();
        ADD_FAILURE() << "expected the tabulated density to be refused";
    } catch (const amrex::RuntimeError& err) {
        EXPECT_NE(
            std::string(err.what())
                .find("density cannot be read from a profile file"),
            std::string::npos)
            << err.what();
    }
}

// A PDE records the fill interpolation of its field, so that an inflow UDF
// registered later uses it too (TKE.interpolation = PiecewiseConstant)
TEST_F(TabulatedProfileTest, a_pde_field_records_its_fill_interpolation)
{
    populate_parameters();
    for (const auto& face : {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    initialize_mesh();

    auto fields = kynema_sgf::pde::create_fields_instance<
        kynema_sgf::pde::TKE, kynema_sgf::fvm::Godunov>(
        time(), mesh().field_repo(),
        kynema_sgf::FieldInterpolator::PiecewiseConstant);
    EXPECT_EQ(
        fields.field.fillpatch_interpolator(),
        kynema_sgf::FieldInterpolator::PiecewiseConstant);
}

TEST_F(TabulatedProfileTest, a_missing_w_column_is_zero)
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

// w on its z faces sits on the z nodes, so the profile is read at k dz
TEST_F(TabulatedProfileTest, w_on_z_faces_is_read_at_the_nodes)
{
    populate_parameters();
    write_profile(
        "tp_w.txt",
        "# z u v w T\n"
        "0.0  0.0  3.0  0.0  300.0\n"
        "8.0 16.0 -5.0  8.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_w.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    // u = 2z, v = 3 - z and w = z, all read at the node heights k dz
    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 1.0_rt},
        {0.0_rt, 3.0_rt, 0.0_rt}, 0.0_rt, 2);
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

// The fill operator of a face-centered field hands its face direction to the
// profile, so that the fill above reads w at the nodes
TEST_F(TabulatedProfileTest, face_fills_pass_their_direction_to_the_profile)
{
    populate_parameters();
    write_profile(
        "tp_w2.txt",
        "# z u v w T\n"
        "0.0  0.0  3.0  0.0  300.0\n"
        "8.0 16.0 -5.0  8.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_w2.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    const kynema_sgf::BCOpCreator<
        kynema_sgf::udf::TabulatedProfile, kynema_sgf::ConstDirichlet>
        creator(vel);
    EXPECT_EQ(creator(2).m_inflow_op.face_dir, 2);
    EXPECT_EQ(creator().m_inflow_op.face_dir, -1);
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
    expect_abort_with(vel, "has 3 columns. Without a header");
}

// The tke boundary condition needs a tke column wherever tke is profiled
TEST_F(TabulatedProfileTest, tke_without_a_tke_column_is_rejected)
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
    initialize_mesh();

    auto& tke = inflow_field("tke", 1, {m_xlo});
    expect_abort_with(tke, "tp_no_tke.txt has no tke column");
}

// Velocity may be profiled from a file without tke: only a tke field filled
// from the profile needs that column
TEST_F(TabulatedProfileTest, velocity_needs_no_tke_column)
{
    populate_parameters();
    write_profile(
        "tp_no_tke2.txt",
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_no_tke2.txt"));
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// interp::linear takes heights closer than the single-precision epsilon for
// one, so a profile with such knots would not be interpolated between them
TEST_F(TabulatedProfileTest, heights_too_close_to_interpolate_are_rejected)
{
    populate_parameters();
    write_profile(
        "tp_close.txt",
        "# z u v T\n"
        "0.0     1.0  0.0  300.0\n"
        "1.0e-8  2.0  0.0  300.0\n"
        "8.0     3.0  0.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_close.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "tp_close.txt line 3 has");
    expect_abort_with(vel, "below which the profile cannot be interpolated");
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
    expect_abort_with(vel, "heights must increase strictly");
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
    expect_abort_with(
        vel,
        "changes sign over the column, so part of the xlo boundary is an "
        "outflow");
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

    {
        // An inflow-outflow face asks for the UDF with its own key
        amrex::ParmParse ppx("xlo");
        ppx.add(
            "velocity.inflow_outflow_type", std::string("TabulatedProfile"));
    }

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

TEST_F(TabulatedProfileTest, outflow_part_of_a_scalar_face_is_extrapolated)
{
    populate_parameters();
    // A transported scalar other than temperature or tke, profiled on an
    // inflow-outflow face
    write_profile(
        "tp_tracer.txt",
        "# z u v T tracer\n"
        "0.0  1.0  0.0  300.0  0.0\n"
        "8.0  1.0  0.0  308.0  16.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_tracer.txt"));
    }
    for (const auto& face : {"ylo", "zlo", "xhi", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("type", std::string("mass_inflow_outflow"));
        pp.add("tracer.inflow_outflow_type", std::string("TabulatedProfile"));
    }
    initialize_mesh();

    // Flow enters through the lower half of xlo and leaves through the upper
    const int kmid = 4;
    const amrex::Real interior = 100.0_rt;
    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    set_inflow_outflow_velocity(vel, kmid);

    auto& tracer = frepo.declare_field("tracer", 1, 1, 1);
    tracer.setVal(interior);
    kynema_sgf::BCScalar bc(tracer);
    // Marked as transported, as the BC setup of its transport equation does
    bc.set_transported();
    bc(0.0_rt);
    kynema_sgf::scalar_bc::register_scalar_dirichlet(
        tracer, mesh(), time(), bc.get_dirichlet_udfs());

    tracer.fillphysbc(0.0_rt);
    tracer.apply_bc_funcs(kynema_sgf::FieldState::New);

    // The profile is held where the flow enters and the interior value is
    // extrapolated where it leaves
    const auto err = xlo_ghost_error(tracer, kmid, 2.0_rt, interior);
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    EXPECT_NEAR(err, 0.0_rt, tol);
}

// On a mass_inflow_outflow face a transported scalar takes the interior value
// where the flow leaves, whatever sets its inflow value. Here the inflow value
// is a constant: the ghost cells start at that constant (0) and the interior
// is 100, so only the outflow half of the face may change
TEST_F(
    TabulatedProfileTest,
    outflow_of_a_constant_transported_scalar_is_extrapolated)
{
    populate_parameters();
    for (const auto& face : {"ylo", "zlo", "xhi", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("type", std::string("mass_inflow_outflow"));
        // The constant inflow value
        pp.add("tracer", 0.0_rt);
    }
    initialize_mesh();

    // Flow enters through the lower half of xlo and leaves through the upper
    const int kmid = 4;
    const amrex::Real interior = 100.0_rt;
    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    set_inflow_outflow_velocity(vel, kmid);

    auto& tracer = frepo.declare_field("tracer", 1, 1, 1);
    tracer.setVal(0.0_rt);
    tracer(0).setVal(interior, 0, 1, 0);
    kynema_sgf::BCScalar bc(tracer);
    bc.set_transported();
    bc(0.0_rt);
    tracer.apply_bc_funcs(kynema_sgf::FieldState::New);

    // 0 (slope 0) where the flow enters, the interior value where it leaves
    const auto err = xlo_ghost_error(tracer, kmid, 0.0_rt, interior);
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    EXPECT_NEAR(err, 0.0_rt, tol);
}

// The BC setup of a transport equation marks its field as transported (the
// production path, not set_transported() called by hand): a passive scalar on
// a mass_inflow_outflow face is extrapolated where the flow leaves, while the
// constant density of ICNS, which is not transported, keeps its value
TEST_F(TabulatedProfileTest, a_transport_equation_marks_its_field)
{
    populate_parameters();
    for (const auto& face : {"ylo", "zlo", "xhi", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("type", std::string("mass_inflow_outflow"));
        amrex::Vector<amrex::Real> uvw{{0.0_rt, 0.0_rt, 0.0_rt}};
        pp.addarr("velocity", uvw);
        pp.add("density", 0.0_rt);
        pp.add("passive_scalar", 0.0_rt);
    }
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("PassiveScalar");

    // Flow enters through the lower half of xlo and leaves through the upper
    const int kmid = 4;
    const amrex::Real interior = 100.0_rt;
    set_inflow_outflow_velocity(sim().repo().get_field("velocity"), kmid);
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Ghost cells start at the inflow constant 0 and the interior is 100
    auto& scalar = sim().repo().get_field("passive_scalar");
    scalar.setVal(0.0_rt);
    scalar(0).setVal(interior, 0, 1, 0);
    scalar.apply_bc_funcs(kynema_sgf::FieldState::New);
    EXPECT_NEAR(xlo_ghost_error(scalar, kmid, 0.0_rt, interior), 0.0_rt, tol);

    auto& density = sim().repo().get_field("density");
    density.setVal(0.0_rt);
    density(0).setVal(interior, 0, 1, 0);
    density.apply_bc_funcs(kynema_sgf::FieldState::New);
    EXPECT_NEAR(xlo_ghost_error(density, 0, 0.0_rt, 0.0_rt), 0.0_rt, tol);
}

// A field that is not transported, such as pressure or a source term, keeps
// its inflow value on the whole face, where the flow leaves too
TEST_F(TabulatedProfileTest, outflow_of_a_field_not_transported_is_kept)
{
    populate_parameters();
    for (const auto& face : {"ylo", "zlo", "xhi", "yhi", "zhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("slip_wall"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("type", std::string("mass_inflow_outflow"));
        // The constant inflow value
        pp.add("other", 0.0_rt);
    }
    initialize_mesh();

    const int kmid = 4;
    const amrex::Real interior = 100.0_rt;
    auto& frepo = mesh().field_repo();
    auto& vel = frepo.declare_field("velocity", 3, 1, 1);
    set_inflow_outflow_velocity(vel, kmid);

    auto& other = frepo.declare_field("other", 1, 1, 1);
    other.setVal(0.0_rt);
    other(0).setVal(interior, 0, 1, 0);
    kynema_sgf::BCScalar bc(other);
    bc(0.0_rt);
    other.apply_bc_funcs(kynema_sgf::FieldState::New);

    // kmid = 0: every ghost cell of the face is expected at 0
    const auto err = xlo_ghost_error(other, 0, 0.0_rt, 0.0_rt);
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    EXPECT_NEAR(err, 0.0_rt, tol);
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
    expect_abort_with(vel, "is directed out of the domain everywhere on xlo");
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

TEST_F(TabulatedProfileTest, only_the_profile_inside_the_domain_is_checked)
{
    populate_parameters();
    // u only turns around at z = 40, so over this 8 m tall domain it falls
    // from 4 to 3.2 and enters everywhere, although the knot above the domain
    // points out of it
    write_profile(
        "tp_far_knot.txt",
        "# z u v T\n"
        "0.0    4.0  0.0  300.0\n"
        "80.0  -4.0  0.0  380.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_far_knot.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

TEST_F(TabulatedProfileTest, a_reversal_between_two_knots_is_found)
{
    populate_parameters();
    // Neither knot lies strictly inside the domain, but u changes sign at
    // z = 5, so only evaluating the ends of the domain catches it
    write_profile(
        "tp_between_knots.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "10.0 -4.0  0.0  310.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_between_knots.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "changes sign over the column, so part of the xlo boundary is an "
        "outflow");
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
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"TerrainDrag"});
    }
    initialize_mesh();

    // The ground is at 3 but no offset was given
    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the ground on xlo is at");
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
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"TerrainDrag"});
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
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"TerrainDrag"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the ground along xlo varies between");
}

// The ground checks concern profiles only: a constant ylo face may stand on
// ground that rises along it while the profiled xlo face stands level
TEST_F(TabulatedProfileTest, a_constant_face_may_stand_on_varying_ground)
{
    populate_parameters();
    write_profile(
        "tp_g6.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    write_terrain_along_x("tp_terrain_xslope.amrwind", 0.0_rt, 4.0_rt);
    {
        amrex::ParmParse pp("xlo");
        pp.add("tabulated_profile_file", std::string("tp_g6.txt"));
    }
    {
        amrex::ParmParse pp("ylo");
        amrex::Vector<amrex::Real> uvw{{7.0_rt, 8.0_rt, 9.0_rt}};
        pp.addarr("velocity", uvw);
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_terrain_xslope.amrwind"));
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"TerrainDrag"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo, m_ylo});
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    const auto err = max_error(
        vel, mesh().Geom(0), profile, m_ylo, {0.0_rt, 0.0_rt, 0.0_rt},
        {7.0_rt, 8.0_rt, 9.0_rt});
    EXPECT_NEAR(err, 0.0_rt, m_tol);
}

// With TerrainDrag active and no TerrainDrag.terrain_file, the ground is the
// TerrainDrag default terrain.amrwind, and it is checked like any other file
TEST_F(TabulatedProfileTest, the_default_terrain_drag_file_is_checked)
{
    populate_parameters();
    write_profile(
        "tp_g8.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    // Rising across the span, so the xlo face does not stand on level ground
    write_terrain("terrain.amrwind", 0.0_rt, 8.0_rt);
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g8.txt"));
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"ABL", "TerrainDrag"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the ground along xlo varies between");
    // The default name is shared with other tests; leave no file behind
    std::remove("terrain.amrwind");
}

namespace {
//! Set up a profiled xlo face over ground rising along it, with the given
//! physics and the terrain in TerrainDrag.terrain_file
void sloped_ground_case(
    const std::string& tag, const amrex::Vector<std::string>& physics)
{
    write_profile(
        "tp_" + tag + ".txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    write_terrain("tp_" + tag + ".amrwind", 0.0_rt, 8.0_rt);
    amrex::ParmParse("TabulatedProfile").add("filename", "tp_" + tag + ".txt");
    amrex::ParmParse("TerrainDrag")
        .add("terrain_file", "tp_" + tag + ".amrwind");
    if (!physics.empty()) {
        amrex::ParmParse("incflo").addarr("physics", physics);
    }
}
} // namespace

// Physics are built in the order listed, and TerrainDrag takes its terrain
// from the waves only when OceanWaves is built before it. Listed after, the
// waves leave TerrainDrag reading its file, which is then checked
TEST_F(TabulatedProfileTest, terrain_drag_before_the_waves_reads_its_file)
{
    populate_parameters();
    sloped_ground_case("wv1", {"ABL", "TerrainDrag", "OceanWaves"});
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the ground along xlo varies between");
}

// OceanWaves built before TerrainDrag gives it wave terrain, and the file is
// not used
TEST_F(TabulatedProfileTest, waves_before_terrain_drag_leave_no_file)
{
    populate_parameters();
    sloped_ground_case("wv2", {"ABL", "OceanWaves", "TerrainDrag"});
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// A vof field built before TerrainDrag keeps it on its file, even with the
// waves; built after, it does not
TEST_F(TabulatedProfileTest, only_a_vof_built_before_terrain_drag_counts)
{
    populate_parameters();
    sloped_ground_case("wv3", {"MultiPhase", "OceanWaves", "TerrainDrag"});
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the ground along xlo varies between");

    amrex::ParmParse("incflo").addarr(
        "physics",
        amrex::Vector<std::string>{"OceanWaves", "TerrainDrag", "MultiPhase"});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});

    // A level set declares no vof field
    amrex::ParmParse("MultiPhase")
        .add("interface_capturing_method", std::string("levelset"));
    amrex::ParmParse("incflo").addarr(
        "physics",
        amrex::Vector<std::string>{"MultiPhase", "OceanWaves", "TerrainDrag"});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// Without TerrainDrag the file is read only by the ABL physics, for a
// terrain-aligned initial profile, and only then is it checked
TEST_F(TabulatedProfileTest, a_terrain_file_is_checked_only_when_read)
{
    populate_parameters();
    sloped_ground_case("wv4", {});
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_wind_profile", true);
        pp.add("terrain_aligned_profile", true);
    }
    initialize_mesh();

    // Neither TerrainDrag nor ABL is built: the file is not used
    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});

    amrex::ParmParse("incflo").addarr(
        "physics", amrex::Vector<std::string>{"ABL"});
    expect_abort_with(vel, "the ground along xlo varies between");
}

// The ABL inputs set the interior only through the ABL physics, so without
// it an offset cannot conflict with them
TEST_F(TabulatedProfileTest, abl_inputs_without_the_abl_physics_are_unused)
{
    populate_parameters();
    write_profile(
        "tp_g13.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g13.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_wind_profile", true);
        pp.add("terrain_aligned_profile", false);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// The ABL physics reads the default terrain.amrwind for a terrain-aligned
// profile when no file is named, so that file is checked, and a missing one
// is an error
TEST_F(TabulatedProfileTest, the_default_abl_terrain_file_is_checked)
{
    populate_parameters();
    write_profile(
        "tp_g14.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    amrex::ParmParse("TabulatedProfile")
        .add("filename", std::string("tp_g14.txt"));
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_wind_profile", true);
        pp.add("terrain_aligned_profile", true);
    }
    amrex::ParmParse("incflo").addarr(
        "physics", amrex::Vector<std::string>{"ABL"});
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    std::remove("terrain.amrwind");
    expect_abort_with(
        vel,
        "cannot open the ABL terrain-aligned profile terrain file "
        "terrain.amrwind");

    write_terrain("terrain.amrwind", 0.0_rt, 8.0_rt);
    expect_abort_with(vel, "the ground along xlo varies between");
    // The default name is shared with other tests; leave no file behind
    std::remove("terrain.amrwind");
}

// With TerrainDrag active, a terrain file that cannot be read is an error
// rather than a check quietly skipped
TEST_F(TabulatedProfileTest, a_missing_terrain_drag_file_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_g9.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g9.txt"));
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_no_such_terrain.amrwind"));
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"ABL", "TerrainDrag"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "cannot open the TerrainDrag terrain file "
        "tp_no_such_terrain.amrwind");
}

TEST_F(TabulatedProfileTest, a_negative_ground_tolerance_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_g10.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g10.txt"));
        pp.add("ground_tolerance", -1.0_rt);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "ground_tolerance must not be negative");
}

// One fill operator serves every inflow face of a field, so another UDF on
// one face cannot be combined with the profile on another
TEST_F(TabulatedProfileTest, another_inflow_udf_on_a_face_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_mix.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_mix.txt"));
    }
    {
        amrex::ParmParse pp("xlo");
        pp.add("velocity.inflow_type", std::string("PowerLawProfile"));
    }
    initialize_mesh();

    // ylo asks for the profile; xlo, also an inflow face, for a power law
    auto& vel = inflow_field("velocity", 3, {m_ylo});
    vel.bc_type()[m_xlo] = BC::mass_inflow;
    expect_abort_with(
        vel,
        "xlo.velocity.inflow_type = PowerLawProfile cannot be combined "
        "with TabulatedProfile");
}

// A top face sits at one height: w changing sign lower down does not make it
// an outflow as long as the flow enters at the top
TEST_F(TabulatedProfileTest, a_top_face_is_checked_at_its_own_height)
{
    populate_parameters();
    write_profile(
        "tp_top.txt",
        "# z u v w T\n"
        "0.0   0.0  0.0   1.0  300.0\n"
        "8.0   0.0  0.0  -1.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_top.txt"));
    }
    initialize_mesh();

    const amrex::Orientation zhi{2, amrex::Orientation::high};
    auto& vel = inflow_field("velocity", 3, {zhi});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// A face uses the profile only when it selects the UDF with the key of its own
// boundary type; other inflow faces take their constant even though
// TabulatedProfile.filename names a profile for every face
TEST_F(
    TabulatedProfileTest, a_face_that_does_not_select_the_profile_is_constant)
{
    populate_parameters();
    write_profile(
        "tp_sel.txt",
        "# z u v T\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_sel.txt"));
    }
    {
        // ylo is a mass_inflow face without velocity.inflow_type
        amrex::ParmParse pp("ylo");
        amrex::Vector<amrex::Real> uvw{{7.0_rt, 8.0_rt, 9.0_rt}};
        pp.addarr("velocity", uvw);
    }
    {
        // yhi is a mass_inflow_outflow face that sets the mass_inflow key,
        // which does not apply to it
        amrex::ParmParse pp("yhi");
        pp.add("velocity.inflow_type", std::string("TabulatedProfile"));
        amrex::Vector<amrex::Real> uvw{{1.0_rt, 2.0_rt, 3.0_rt}};
        pp.addarr("velocity", uvw);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    vel.bc_type()[m_ylo] = BC::mass_inflow;
    const amrex::Orientation yhi{1, amrex::Orientation::high};
    vel.bc_type()[yhi] = BC::mass_inflow_outflow;
    const kynema_sgf::udf::TabulatedProfile profile(vel);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> zero{
        0.0_rt, 0.0_rt, 0.0_rt};
    EXPECT_NEAR(
        max_error(
            vel, mesh().Geom(0), profile, m_xlo, {2.0_rt, -1.0_rt, 0.0_rt},
            {0.0_rt, 3.0_rt, 0.0_rt}),
        0.0_rt, m_tol);
    EXPECT_NEAR(
        max_error(
            vel, mesh().Geom(0), profile, m_ylo, zero,
            {7.0_rt, 8.0_rt, 9.0_rt}),
        0.0_rt, m_tol);
    EXPECT_NEAR(
        max_error(
            vel, mesh().Geom(0), profile, yhi, zero, {1.0_rt, 2.0_rt, 3.0_rt}),
        0.0_rt, m_tol);
}

// A scalar face that does not use the profile needs its constant: there is
// no sensible default (0 K for temperature)
TEST_F(TabulatedProfileTest, a_scalar_face_without_a_constant_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_sc.txt",
        "# z u v T\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_sc.txt"));
    }
    initialize_mesh();

    auto& temp = inflow_field("temperature", 1, {m_xlo});
    temp.bc_type()[m_ylo] = BC::mass_inflow;
    expect_abort_with(
        temp,
        "ylo does not use the profile for temperature, so it needs the "
        "constant value ylo.temperature");
}

// A raised ground tolerance accepts a face whose ground varies within it
TEST_F(TabulatedProfileTest, the_ground_tolerance_can_be_raised)
{
    populate_parameters();
    write_profile(
        "tp_g11.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    // The xlo ground rises from 0 to 1.5 along y, more than the default
    // tolerance (dz = 1); its mean is 0.75
    write_terrain("tp_terrain_gentle.amrwind", 0.0_rt, 1.5_rt);
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g11.txt"));
        pp.add("zoffset", 0.75_rt);
        pp.add("ground_tolerance", 2.0_rt);
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_terrain_gentle.amrwind"));
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"TerrainDrag"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// With the initial profile aligned with the terrain, an offset is consistent
TEST_F(
    TabulatedProfileTest, an_offset_with_a_terrain_aligned_interior_is_accepted)
{
    populate_parameters();
    write_profile(
        "tp_g12.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g12.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    // The ABL physics reads the terrain for the aligned profile, so the
    // offset is checked against it
    write_terrain("tp_terrain_flat3.amrwind", 3.0_rt, 3.0_rt);
    amrex::ParmParse("TerrainDrag")
        .add("terrain_file", std::string("tp_terrain_flat3.amrwind"));
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_wind_profile", true);
        pp.add("terrain_aligned_profile", true);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"ABL"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// On raised ground the part of the column below the ground is buried in the
// terrain: a profile pointing out of the domain only there is accepted
TEST_F(TabulatedProfileTest, the_buried_part_of_a_raised_face_is_not_checked)
{
    populate_parameters();
    write_profile(
        "tp_buried.txt",
        "# z u v T\n"
        "-3.0  -1.0  0.0  300.0\n"
        "-1.0  -1.0  0.0  300.0\n"
        "0.0    1.0  0.0  300.0\n"
        "5.0    1.0  0.0  305.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_buried.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
}

// Above the ground the check still applies
TEST_F(TabulatedProfileTest, a_raised_face_reversing_above_ground_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_buried2.txt",
        "# z u v T\n"
        "0.0    1.0  0.0  300.0\n"
        "5.0   -1.0  0.0  305.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_buried2.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "changes sign over the column");
}

// A top face the flow leaves through at the top is refused as mass_inflow
TEST_F(TabulatedProfileTest, a_top_face_leaving_at_the_top_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_top2.txt",
        "# z u v w T\n"
        "0.0   0.0  0.0  -1.0  300.0\n"
        "8.0   0.0  0.0   1.0  308.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_top2.txt"));
    }
    initialize_mesh();

    const amrex::Orientation zhi{2, amrex::Orientation::high};
    auto& vel = inflow_field("velocity", 3, {zhi});
    expect_abort_with(vel, "is directed out of the domain everywhere on zhi");
}

// Only velocity and single-component scalars can be read from a profile
TEST_F(TabulatedProfileTest, a_multicomponent_scalar_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_mc.txt",
        "# z u v T\n"
        "0.0  0.0  3.0  300.0\n"
        "8.0 16.0 -5.0  308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_mc.txt"));
    initialize_mesh();

    auto& pair = inflow_field("pair", 2, {m_xlo});
    expect_abort_with(
        pair, "only velocity and single-component scalars can be read");
}

// Density is not tabulated, even when its column is in the file
TEST_F(TabulatedProfileTest, a_tabulated_density_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_rho.txt",
        "# z u v T density\n"
        "0.0  0.0  3.0  300.0  1.2\n"
        "8.0 16.0 -5.0  308.0  1.1\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_rho.txt"));
    initialize_mesh();

    auto& rho = inflow_field("density", 1, {m_xlo});
    expect_abort_with(rho, "density cannot be read from a profile file");
}

// The UDF was requested but no face both selects it and has a file
TEST_F(TabulatedProfileTest, a_profile_with_no_file_is_rejected)
{
    populate_parameters();
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "but no inflow face uses it with a profile file");
}

// A bottom or top face holds the value at its own height: a ghost cell of a
// cell-centered field and the face node of w both take the profile there
TEST_F(TabulatedProfileTest, bottom_and_top_faces_are_read_at_their_height)
{
    populate_parameters();
    write_profile(
        "tp_ztop.txt",
        "# z u v w T\n"
        "-8.0  0.0  0.0  0.0  292.0\n"
        "16.0  0.0  0.0  3.0  316.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_ztop.txt"));
    }
    initialize_mesh();

    const amrex::Orientation zlo{2, amrex::Orientation::low};
    const amrex::Orientation zhi{2, amrex::Orientation::high};
    const auto& geom = mesh().Geom(0);
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // T = 300 + z, tabulated beyond the domain: 300 at the bottom face and
    // 308 at the top face (the domain is 8 high), not 299.5 and 308.5 at the
    // ghost cell centers
    auto& temp = inflow_field("temperature", 1, {zlo, zhi});
    const kynema_sgf::udf::TabulatedProfile tprof(temp);
    EXPECT_NEAR(
        value_at(tprof, geom, zlo, amrex::IntVect{2, 2, -1}, 0), 300.0_rt,
        tol * 300.0_rt);
    EXPECT_NEAR(
        value_at(tprof, geom, zhi, amrex::IntVect{2, 2, 8}, 0), 308.0_rt,
        tol * 300.0_rt);

    // w = 1 + z/8 as set on the bottom face (index 0, the fill of the MAC
    // inflow) is 1, not 1.0625 half a cell up
    auto& vel = inflow_field("velocity", 3, {zlo});
    const kynema_sgf::udf::TabulatedProfile vprof(vel);
    EXPECT_NEAR(
        value_at(vprof, geom, zlo, amrex::IntVect{2, 2, 0}, 2), 1.0_rt, tol);
}

// A bump between two cell centers (3.5 and 4.5) is still ground that varies
TEST_F(TabulatedProfileTest, ground_varying_between_cell_centers_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_g5.txt",
        "# z u v T\n"
        "0.0   4.0  0.0  300.0\n"
        "8.0   4.0  0.0  308.0\n");
    write_terrain_along_y(
        "tp_terrain_bump.amrwind", {0.0_rt, 3.5_rt, 4.0_rt, 4.5_rt, 8.0_rt},
        {0.0_rt, 0.0_rt, 8.0_rt, 0.0_rt, 0.0_rt});
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g5.txt"));
    }
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", std::string("tp_terrain_bump.amrwind"));
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"TerrainDrag"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the ground along xlo varies between");
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
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"ABL"});
    }
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "ABL.terrain_aligned_profile is off");
}

// ABL.initial_wind_profile only sets velocity, temperature and tke, so the
// offset of another scalar cannot conflict with it
TEST_F(TabulatedProfileTest, a_scalar_offset_ignores_the_abl_initial_profile)
{
    populate_parameters();
    write_profile(
        "tp_g7.txt",
        "# z u v T tracer\n"
        "0.0   4.0  0.0  300.0  0.0\n"
        "8.0   4.0  0.0  308.0  1.0\n");
    {
        amrex::ParmParse pp("TabulatedProfile");
        pp.add("filename", std::string("tp_g7.txt"));
        pp.add("zoffset", 3.0_rt);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_wind_profile", true);
        pp.add("terrain_aligned_profile", false);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.addarr("physics", amrex::Vector<std::string>{"ABL"});
    }
    initialize_mesh();

    auto& tracer = inflow_field("tracer", 1, {m_xlo});
    EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{tracer});
}

TEST_F(TabulatedProfileTest, a_trailing_word_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e1.txt",
        "# z u v T\n"
        "0.0 1.0 2.0 300.0 junk\n"
        "8.0 3.0 4.0 308.0 junk\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e1.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel, "tp_e1.txt line 2, column 5: 'junk' is not a number");
}

TEST_F(TabulatedProfileTest, a_word_where_a_number_belongs_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e2.txt",
        "# z u v T\n"
        "0.0 abc 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e2.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "tp_e2.txt line 2, column 2: 'abc' is not a number");
}

TEST_F(TabulatedProfileTest, a_nan_in_the_file_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e3.txt",
        "# z u v T\n"
        "0.0 nan 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e3.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel, "tp_e3.txt line 2, column 2: 'nan' is not a finite value");
}

TEST_F(TabulatedProfileTest, an_infinity_in_the_file_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e4.txt",
        "# z u v T\n"
        "0.0 inf 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e4.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel, "tp_e4.txt line 2, column 2: 'inf' is not a finite value");
}

TEST_F(TabulatedProfileTest, an_unrepresentable_number_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e5.txt",
        "# z u v T\n"
        "0.0 1e400 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e5.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel,
        "tp_e5.txt line 2, column 2: '1e400' is too large or too small to "
        "represent");
}

TEST_F(TabulatedProfileTest, rows_of_different_widths_are_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e6.txt",
        "# z u v T\n"
        "0.0 1.0 2.0 300.0\n"
        "8.0 3.0 4.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e6.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(
        vel, "tp_e6.txt line 3 has 3 columns but tp_e6.txt line 2 has 4");
}

TEST_F(TabulatedProfileTest, a_file_with_no_data_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e7.txt",
        "# z u v T\n"
        "\n"
        "# nothing follows\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e7.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "tp_e7.txt holds no profile data");
}

TEST_F(TabulatedProfileTest, a_single_height_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e8.txt",
        "# z u v T\n"
        "0.0 1.0 2.0 300.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e8.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "tp_e8.txt must tabulate at least two heights");
}

TEST_F(TabulatedProfileTest, a_repeated_column_name_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e9.txt",
        "# z u u T\n"
        "0.0 1.0 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e9.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the header of tp_e9.txt names 'u' more than once");
}

TEST_F(TabulatedProfileTest, a_height_with_no_values_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e10.txt",
        "# z u v T\n"
        "0.0\n"
        "8.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e10.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "tp_e10.txt line 2 holds a height and no values");
}

TEST_F(TabulatedProfileTest, a_repeated_height_name_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e11.txt",
        "# z z u v T\n"
        "0.0 0.0 1.0 2.0 300.0\n"
        "8.0 8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e11.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "the header of tp_e11.txt names 'z' more than once");
}

// 1e40 and 1e-40 are finite doubles but overflow and underflow a float
TEST_F(TabulatedProfileTest, a_number_too_large_for_real_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e12.txt",
        "# z u v T\n"
        "0.0 1e40 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e12.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    if (std::numeric_limits<amrex::Real>::max() < 1.0e40) {
        expect_abort_with(
            vel,
            "tp_e12.txt line 2, column 2: '1e40' is too large or too small to "
            "represent");
    } else {
        EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
    }
}

TEST_F(TabulatedProfileTest, a_number_too_small_for_real_is_rejected)
{
    populate_parameters();
    write_profile(
        "tp_e13.txt",
        "# z u v T\n"
        "0.0 1e-40 2.0 300.0\n"
        "8.0 3.0 4.0 308.0\n");
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_e13.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    if (std::numeric_limits<amrex::Real>::min() > 1.0e-40) {
        expect_abort_with(
            vel,
            "tp_e13.txt line 2, column 2: '1e-40' is too large or too small "
            "to represent");
    } else {
        EXPECT_NO_THROW(kynema_sgf::udf::TabulatedProfile{vel});
    }
}

// The same file under two spellings of its name is still the RANS file, and
// the two inputs cannot share it
TEST_F(TabulatedProfileTest, the_rans_profile_file_cannot_be_shared)
{
    populate_parameters();
    write_profile(
        "tp_rans.txt",
        "0.0 1.0 2.0 0.0 0.5\n"
        "8.0 3.0 4.0 0.0 0.4\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("rans_1dprofile_file", std::string("./tp_rans.txt"));
    }
    amrex::ParmParse pp("TabulatedProfile");
    pp.add("filename", std::string("tp_rans.txt"));
    initialize_mesh();

    auto& vel = inflow_field("velocity", 3, {m_xlo});
    expect_abort_with(vel, "is also used as ABL.rans_1dprofile_file");
}

} // namespace kynema_sgf_tests
