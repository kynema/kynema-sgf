#include "src/physics/udfs/TabulatedProfile.H"
#include "src/core/Field.H"
#include "src/core/FieldRepo.H"
#include "src/incflo_enums.H"
#include "src/utilities/constants.H"
#include "src/utilities/io_utils.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_REAL.H"
#include "AMReX_Utility.H"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <system_error>
#include <utility>

using namespace amrex::literals;

namespace kynema_sgf::udf {

namespace {

//! Face names in the order given by amrex::Orientation
const amrex::Vector<std::string> face_names = {"xlo", "ylo", "zlo",
                                               "xhi", "yhi", "zhi"};

//! A vertical profile read from one file
struct ProfileData
{
    //! Tabulated heights, strictly increasing
    amrex::Vector<amrex::Real> z;

    //! Names of the columns after the height column
    amrex::Vector<std::string> colnames;

    //! Column values, each of the same length as z
    amrex::Vector<amrex::Vector<amrex::Real>> cols;

    //! Whether the column names came from a header rather than being assumed
    bool has_header{false};

    //! Line of the file the header was taken from, when there is one
    int header_line{0};

    [[nodiscard]] int column_index(const std::string& name) const
    {
        for (int i = 0; i < static_cast<int>(colnames.size()); ++i) {
            if (colnames[i] == name) {
                return i;
            }
        }
        return -1;
    }
};

/** Whether two file names point at the same file
 *
 *  Compares the files themselves, so './a.txt', 'a.txt' and a link to it are
 *  one file. Falls back to the normalized paths when a file cannot be found.
 *
 *  \param a First file name
 *  \param b Second file name
 *  \return True when both names refer to the same file
 */
bool same_file(const std::string& a, const std::string& b)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    const bool equal = fs::equivalent(a, b, ec);
    if (!ec) {
        return equal;
    }
    const auto pa = fs::weakly_canonical(a, ec);
    const auto pb = ec ? fs::path() : fs::weakly_canonical(b, ec);
    if (!ec) {
        return pa == pb;
    }
    return fs::path(a).lexically_normal() == fs::path(b).lexically_normal();
}

/** Map the spellings accepted in a header onto the field names used internally
 *
 *  \param name Name as written
 *  \return Lower-case name, with T, theta and temp mapped to temperature
 */
std::string canonical_name(const std::string& name)
{
    auto lname = amrex::toLower(name);
    if ((lname == "t") || (lname == "theta") || (lname == "temp")) {
        return "temperature";
    }
    return lname;
}

/** Assume column names for a file that carries no header
 *
 *  Four columns are ``z u v T`` and five are ``z u v T tke``; any other width
 *  is ambiguous and must be labeled by a header instead.
 *
 *  \param ncols Number of columns of the data, height included
 *  \param fname Profile file, for the message
 *  \return Names of the columns after the height
 */
amrex::Vector<std::string>
assumed_columns(const int ncols, const std::string& fname)
{
    if (ncols == 4) {
        return {"u", "v", "temperature"};
    }
    if (ncols == 5) {
        return {"u", "v", "temperature", "tke"};
    }
    amrex::Abort(
        "TabulatedProfile: " + fname + " has " + std::to_string(ncols) +
        " columns. Without a header only 4 columns (z u v T) or 5 columns "
        "(z u v T tke) can be interpreted. Add a header line, for example "
        "'# z u v T tke'.");
    return {};
}

/** The words of a comment line, as column names
 *
 *  Leading '#' characters are dropped and commas count as spaces, so that
 *  '## z u v T' and '# z, u, v, T' read like '# z u v T'.
 *
 *  \param line Comment line, starting with '#'
 *  \return Its words, in canonical form
 */
amrex::Vector<std::string> comment_names(std::string line)
{
    std::ranges::replace(line, ',', ' ');
    const auto start = line.find_first_not_of("# \t");
    std::istringstream iss(
        (start == std::string::npos) ? std::string() : line.substr(start));
    amrex::Vector<std::string> names;
    std::string name;
    while (iss >> name) {
        names.push_back(canonical_name(name));
    }
    return names;
}

/** Reject a header that names a column more than once
 *
 *  The height column counts too: '# z z u v T' names z twice.
 *
 *  \param names Column names of the header, height first
 *  \param fname Profile file, for the message
 *  \param lineno Line the header was read from
 */
void check_unique_names(
    const amrex::Vector<std::string>& names,
    const std::string& fname,
    const int lineno)
{
    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(names.size()); ++j) {
            if (names[i] == names[j]) {
                amrex::Abort(
                    "TabulatedProfile: the header of " + fname + " names '" +
                    names[i] + "' more than once (line " +
                    std::to_string(lineno) +
                    "; if that line is a note, reword it so that it does not "
                    "start with 'z')");
            }
        }
    }
}

/** How many column names this reader knows a comment holds
 *
 *  Used to tell a header from a note: a header names u, v, w, temperature or
 *  tke, a note rarely does.
 *
 *  \param names Words of the comment, in canonical form
 *  \return Number of words that are u, v, w, temperature or tke
 */
int known_columns(const amrex::Vector<std::string>& names)
{
    int count = 0;
    for (const auto& n : names) {
        if ((n == "u") || (n == "v") || (n == "w") || (n == "temperature") ||
            (n == "tke")) {
            ++count;
        }
    }
    return count;
}

/** A value written with every digit it holds, so that two close heights do
 *  not print the same in a message
 *
 *  \param val Value to write
 *  \return The value as text
 */
std::string precise(const amrex::Real val)
{
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10)
        << val;
    return oss.str();
}

/** Where a problem was found, for a message the reader can act on
 *
 *  \param fname Profile file
 *  \param lineno Line number
 *  \param col Column number, or 0 for the whole line
 *  \return Text such as "file line 3, column 2"
 */
std::string at(const std::string& fname, const int lineno, const int col = 0)
{
    auto where = fname + " line " + std::to_string(lineno);
    if (col > 0) {
        where += ", column " + std::to_string(col);
    }
    return where;
}

/** Read one number, saying exactly what is wrong when it is not one
 *
 *  Stream extraction stops at the first character it cannot use, which quietly
 *  drops the rest of a line and turns a typo into a puzzling complaint about
 *  the number of columns. Parse the whole token instead.
 *
 *  \param tok Token to read
 *  \param fname Profile file, for the message
 *  \param lineno Line of the token
 *  \param col Column of the token
 *  \return The value
 */
amrex::Real parse_value(
    const std::string& tok,
    const std::string& fname,
    const int lineno,
    const int col)
{
    const char* begin = tok.c_str();
    char* end = nullptr;
    errno = 0;
    const double val = std::strtod(begin, &end);

    if (end != begin + tok.size()) {
        amrex::Abort(
            "TabulatedProfile: " + at(fname, lineno, col) + ": '" + tok +
            "' is not a number");
    }
    if (errno == ERANGE) {
        amrex::Abort(
            "TabulatedProfile: " + at(fname, lineno, col) + ": '" + tok +
            "' is too large or too small to represent");
    }
    if (!std::isfinite(val)) {
        amrex::Abort(
            "TabulatedProfile: " + at(fname, lineno, col) + ": '" + tok +
            "' is not a finite value");
    }
    // strtod checks the range of a double; in single precision a finite
    // double such as 1e40 or 1e-40 still overflows or underflows amrex::Real
    const double mag = std::abs(val);
    constexpr double real_max = std::numeric_limits<amrex::Real>::max();
    constexpr double real_min = std::numeric_limits<amrex::Real>::min();
    if ((mag > real_max) || ((mag > 0.0) && (mag < real_min))) {
        amrex::Abort(
            "TabulatedProfile: " + at(fname, lineno, col) + ": '" + tok +
            "' is too large or too small to represent");
    }
    return static_cast<amrex::Real>(val);
}

/** Read a whitespace-separated profile file
 *
 *  \param fname Profile file
 *  \return Heights, column names and columns of the profile
 */
ProfileData read_profile_file(const std::string& fname)
{
    std::ifstream infile(fname, std::ios::in);
    if (!infile.good()) {
        amrex::Abort("TabulatedProfile: cannot open profile file " + fname);
    }

    ProfileData prof;
    amrex::Vector<amrex::Vector<amrex::Real>> rows;
    amrex::Vector<int> row_lines;
    // Comment lines ahead of the data that could be the header, with the
    // line each came from
    amrex::Vector<amrex::Vector<std::string>> headers;
    amrex::Vector<int> header_lines;
    // Comments ahead of the data that name a known column but do not start
    // with z, such as '# height u v w T', with their lines
    amrex::Vector<amrex::Vector<std::string>> lookalikes;
    amrex::Vector<int> lookalike_lines;
    std::string line;
    int lineno = 0;

    while (std::getline(infile, line)) {
        ++lineno;
        // A UTF-8 byte order mark ahead of the first line is not content
        if ((lineno == 1) && line.starts_with("\xEF\xBB\xBF")) {
            line.erase(0, 3);
        }
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            continue;
        }
        if (line[first] == '#') {
            // Only a header ahead of the data names the columns. Every
            // comment there that could be one is kept, and whether the file
            // has a header at all is decided once the width of the data is
            // known, since a note can also start with z
            if (rows.empty()) {
                const auto names = comment_names(line);
                if (!names.empty() && (names[0] == "z")) {
                    headers.push_back(names);
                    header_lines.push_back(lineno);
                } else if (known_columns(names) > 0) {
                    lookalikes.push_back(names);
                    lookalike_lines.push_back(lineno);
                }
            }
            continue;
        }

        std::istringstream iss(line);
        amrex::Vector<std::string> tokens;
        std::string tok;
        while (iss >> tok) {
            tokens.push_back(tok);
        }
        if (tokens.empty()) {
            continue;
        }

        amrex::Vector<amrex::Real> row;
        row.reserve(tokens.size());
        for (int c = 0; c < static_cast<int>(tokens.size()); ++c) {
            row.push_back(parse_value(tokens[c], fname, lineno, c + 1));
        }

        if (row.size() < 2) {
            amrex::Abort(
                "TabulatedProfile: " + at(fname, lineno) +
                " holds a height and no values");
        }
        if (!rows.empty() && (row.size() != rows[0].size())) {
            amrex::Abort(
                "TabulatedProfile: " + at(fname, lineno) + " has " +
                std::to_string(row.size()) + " columns but " +
                at(fname, row_lines[0]) + " has " +
                std::to_string(rows[0].size()));
        }
        rows.push_back(row);
        row_lines.push_back(lineno);
    }
    infile.close();

    if (rows.empty()) {
        amrex::Abort("TabulatedProfile: " + fname + " holds no profile data");
    }
    if (rows.size() < 2) {
        amrex::Abort(
            "TabulatedProfile: " + fname +
            " must tabulate at least two heights, but holds only one");
    }

    const int ncols = static_cast<int>(rows[0].size());

    // The header is the last candidate as wide as the data. A candidate that
    // is not can be a note ('# z is the height above ground') or a header that
    // does not fit the data. It is taken as a wrong header when it names a
    // column this reader knows, so that a mistake is reported rather than the
    // columns being guessed; otherwise it is a note and the file is read as
    // headerless.
    int header_idx = -1;
    for (int i = 0; i < static_cast<int>(headers.size()); ++i) {
        if (static_cast<int>(headers[i].size()) == ncols) {
            header_idx = i;
        }
    }
    if (header_idx < 0) {
        for (int i = 0; i < static_cast<int>(headers.size()); ++i) {
            if (known_columns(headers[i]) > 0) {
                amrex::Abort(
                    "TabulatedProfile: " + at(fname, header_lines[i]) +
                    " reads as a header naming " +
                    std::to_string(headers[i].size()) +
                    " columns, but the data has " + std::to_string(ncols) +
                    ". Fix the header, or reword the comment so that it does "
                    "not start with 'z'.");
            }
        }
    }

    prof.has_header = (header_idx >= 0);
    if (prof.has_header) {
        const auto& header = headers[header_idx];
        prof.header_line = header_lines[header_idx];
        check_unique_names(header, fname, prof.header_line);
        prof.colnames.assign(header.begin() + 1, header.end());
    } else {
        // A comment as wide as the data that names a known column but does
        // not start with z ('# height u v w T', or '# height u humidity speed'
        // with a single one; '## z u v w T' is fine) would otherwise be
        // skipped and the columns guessed
        for (int i = 0; i < static_cast<int>(lookalikes.size()); ++i) {
            if (static_cast<int>(lookalikes[i].size()) == ncols) {
                amrex::Abort(
                    "TabulatedProfile: " + at(fname, lookalike_lines[i]) +
                    " looks like a header for the " + std::to_string(ncols) +
                    " columns of the data, but a header must start with 'z'");
            }
        }
        prof.colnames = assumed_columns(ncols, fname);
    }

    const int nz = static_cast<int>(rows.size());
    prof.z.resize(nz);
    prof.cols.resize(ncols - 1);
    for (auto& col : prof.cols) {
        col.resize(nz);
    }
    // interp::linear takes two heights closer than this for one and returns
    // the upper value between them, in either precision
    constexpr auto min_spacing =
        static_cast<amrex::Real>(std::numeric_limits<float>::epsilon());
    for (int k = 0; k < nz; ++k) {
        prof.z[k] = rows[k][0];
        if ((k > 0) && (prof.z[k] <= prof.z[k - 1])) {
            amrex::Abort(
                "TabulatedProfile: heights must increase strictly, but " +
                at(fname, row_lines[k]) + " has " + precise(prof.z[k]) +
                " after " + precise(prof.z[k - 1]));
        }
        if ((k > 0) && ((prof.z[k] - prof.z[k - 1]) <= min_spacing)) {
            amrex::Abort(
                "TabulatedProfile: " + at(fname, row_lines[k]) + " has " +
                precise(prof.z[k]) + " after " + precise(prof.z[k - 1]) +
                ", closer than " + precise(min_spacing) +
                ", below which the profile cannot be interpolated");
        }
        for (int c = 0; c < ncols - 1; ++c) {
            prof.cols[c][k] = rows[k][c + 1];
        }
    }

    return prof;
}

/** Names of the columns supplying each component of a field
 *
 *  Velocity draws on ``u``, ``v`` and ``w``; every other field draws on a
 *  column named after the field itself, so that scalars added later are picked
 *  up without touching the reader.
 *
 *  \param field_name Name of the field
 *  \param ncomp Number of components of the field
 *  \return Column name for each component
 */
amrex::Vector<std::string>
wanted_columns(const std::string& field_name, const int ncomp)
{
    if (field_name == "velocity") {
        amrex::Vector<std::string> names = {"u", "v", "w"};
        names.resize(ncomp);
        return names;
    }
    amrex::Vector<std::string> names(ncomp, canonical_name(field_name));
    return names;
}

/** Check that a pure inflow face really does have flow entering everywhere
 *
 *  A veering profile can reverse the normal component partway up the column,
 *  which leaves part of a ``mass_inflow`` face acting as an outflow. That is
 *  what ``mass_inflow_outflow`` is for, so say so rather than injecting flow
 *  backwards through the boundary.
 *
 *  \param heights Heights of every face profile, concatenated
 *  \param vals Values of every face profile, concatenated and interleaved
 *  \param offset Index at which this face profile begins
 *  \param nz Number of heights of this face profile
 *  \param ncomp Number of components of the field
 *  \param face Face index, in amrex::Orientation order
 *  \param zlo Bottom of the domain, measured from the ground of the face
 *  \param zhi Top of the domain, measured from the ground of the face
 *  \param fname Profile file, for the message
 */
void check_inflow_direction(
    const amrex::Vector<amrex::Real>& heights,
    const amrex::Vector<amrex::Real>& vals,
    const int offset,
    const int nz,
    const int ncomp,
    const int face,
    const amrex::Real zlo,
    const amrex::Real zhi,
    const std::string& fname)
{
    const int dir = face % AMREX_SPACEDIM;
    const bool is_low = (face < AMREX_SPACEDIM);

    // Flow enters through a low face when the normal component is positive
    const amrex::Real into = is_low ? 1.0_rt : -1.0_rt;

    // Only the part of the column the domain actually reaches matters. The
    // profile is linear between entries and held constant beyond the table,
    // so over [zlo, zhi] its extremes lie at the two ends of that range or at
    // an entry inside it. Those are evaluated just as the boundary fill does
    const auto* zbegin = heights.data() + offset;
    const auto* zend = zbegin + nz;
    const auto* ybegin =
        vals.data() + (static_cast<std::ptrdiff_t>(ncomp) * offset);

    bool enters = false;
    bool leaves = false;
    const auto classify = [&](const amrex::Real un) {
        if (un > constants::TIGHT_TOL) {
            enters = true;
        }
        if (un < -constants::TIGHT_TOL) {
            leaves = true;
        }
    };
    if (dir == AMREX_SPACEDIM - 1) {
        // A bottom or top face sits at one height, so only the profile there
        // decides whether the flow enters through it
        classify(
            into * interp::linear(
                       zbegin, zend, ybegin, is_low ? zlo : zhi, ncomp, dir));
    } else {
        classify(into * interp::linear(zbegin, zend, ybegin, zlo, ncomp, dir));
        classify(into * interp::linear(zbegin, zend, ybegin, zhi, ncomp, dir));
        for (int k = 0; k < nz; ++k) {
            const auto zk = heights[offset + k];
            if ((zk > zlo) && (zk < zhi)) {
                classify(into * vals[(ncomp * (offset + k)) + dir]);
            }
        }
    }

    if (leaves && enters) {
        amrex::Abort(
            "TabulatedProfile: the normal velocity tabulated in " + fname +
            " changes sign over the column, so part of the " +
            face_names[face] + " boundary is an outflow. Set " +
            face_names[face] +
            ".type = mass_inflow_outflow rather than mass_inflow.");
    }
    if (leaves) {
        amrex::Abort(
            "TabulatedProfile: the normal velocity tabulated in " + fname +
            " is directed out of the domain everywhere on " + face_names[face] +
            ", which is declared mass_inflow.");
    }
}

/** Check the ground along a face against the offset the profile was given
 *
 *  Only a uniform lift is supported, so the ground along an inflow face has to
 *  be flat, and the offset has to be the height it sits at. Ground that varies
 *  along the face would vary the inflow area with it, and the inflow-outflow
 *  solvability correction would then rescale the profile that was asked for.
 *
 *  \param xterrain Terrain grid x coordinates
 *  \param yterrain Terrain grid y coordinates
 *  \param zterrain Terrain heights, indexed [i * ny + j]
 *  \param geom Level-0 geometry
 *  \param face Face index, in amrex::Orientation order
 *  \param zoffset Ground height the profile on this face was given
 *  \param tol Tolerance on the ground variation and on the offset
 */
void check_ground_height(
    const amrex::Vector<amrex::Real>& xterrain,
    const amrex::Vector<amrex::Real>& yterrain,
    const amrex::Vector<amrex::Real>& zterrain,
    const amrex::Geometry& geom,
    const int face,
    const amrex::Real zoffset,
    const amrex::Real tol)
{
    const int dir = face % AMREX_SPACEDIM;
    if (dir == AMREX_SPACEDIM - 1) {
        return;
    }

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();

    // Along the face the bilinear ground is linear between terrain knots, so
    // its extremes are at the ends of the face or at the knots on it. Sampling
    // those, rather than cell centers, does not depend on the mesh or on how
    // the boundary is refined.
    const int tdir = (dir == 0) ? 1 : 0;
    const amrex::Real ncoord =
        (face < AMREX_SPACEDIM) ? problo[dir] : probhi[dir];
    const auto& tknots = (tdir == 0) ? xterrain : yterrain;
    amrex::Vector<amrex::Real> tcoords{problo[tdir], probhi[tdir]};
    for (const auto t : tknots) {
        if ((t > problo[tdir]) && (t < probhi[tdir])) {
            tcoords.push_back(t);
        }
    }

    amrex::Real zmin = constants::LARGE_NUM;
    amrex::Real zmax = -constants::LARGE_NUM;
    for (const auto tcoord : tcoords) {
        const auto xco = (dir == 0) ? ncoord : tcoord;
        const auto yco = (dir == 0) ? tcoord : ncoord;
        const auto zg = interp::bilinear(
            xterrain.data(), xterrain.data() + xterrain.size(), yterrain.data(),
            yterrain.data() + yterrain.size(), zterrain.data(), xco, yco);
        zmin = amrex::min(zmin, zg);
        zmax = amrex::max(zmax, zg);
    }

    if ((zmax - zmin) > tol) {
        amrex::Abort(
            "TabulatedProfile: the ground along " + face_names[face] +
            " varies between " + precise(zmin) + " and " + precise(zmax) +
            ", and only a uniform lift is supported. Raise "
            "TabulatedProfile.ground_tolerance to accept this face, or drive "
            "it with a boundary plane instead.");
    }

    const auto zground = 0.5_rt * (zmin + zmax);
    if (std::abs(zground - zoffset) > tol) {
        amrex::Abort(
            "TabulatedProfile: the ground on " + face_names[face] + " is at " +
            precise(zground) +
            " but the profile is offset "
            "by " +
            precise(zoffset) +
            ". Set the offset to the ground height so that the profile and the "
            "interior are measured from the same place.");
    }
}

} // namespace

TabulatedProfile::TabulatedProfile(const Field& fld)
{
    // This capability is activated with the following in the input file:
    // xlo.type = "mass_inflow"
    // xlo.velocity.inflow_type = TabulatedProfile
    // TabulatedProfile.filename = inflow_profile.txt

    // Density is not tabulated: every inflow face gives it as a constant
    if (fld.name() == "density") {
        amrex::Abort(
            "TabulatedProfile: density cannot be read from a profile file; "
            "give every inflow face a constant <face>.density instead");
    }

    const int ncomp = fld.num_comp();
    // Velocity reads u, v and w; any other field reads the one column named
    // after it, so it must be a single-component scalar
    if ((fld.name() != "velocity") && (ncomp != 1)) {
        amrex::Abort(
            "TabulatedProfile: " + fld.name() + " has " +
            std::to_string(ncomp) +
            " components, but only velocity and single-component scalars can "
            "be read from a profile file");
    }
    AMREX_ALWAYS_ASSERT(ncomp <= AMREX_SPACEDIM);
    m_op.ncomp = ncomp;

    amrex::ParmParse pp("TabulatedProfile");
    std::string default_file;
    pp.query("filename", default_file);

    // Heights in the file are measured from here rather than from the bottom
    // of the domain, which lets a profile given above ground be used on a
    // boundary that sits on uniformly raised ground
    amrex::Real default_zoffset = 0.0_rt;
    pp.query("zoffset", default_zoffset);

    // The checks below follow what the run does with its inputs, so they are
    // decided from the physics it builds, in the order they are listed
    amrex::Vector<std::string> physics;
    amrex::ParmParse("incflo").queryarr("physics", physics);
    const auto has_physics = [&physics](const std::string& name) {
        return std::ranges::find(physics, name) != physics.end();
    };

    // The interior profile and the boundary have to be measured from the same
    // place, or the inflow fights the interior at the boundary. The ABL
    // physics sets the interior from ABL.initial_wind_profile; without it
    // the ABL inputs are not used.
    const bool abl = has_physics("ABL");
    amrex::ParmParse pp_abl("ABL");
    bool init_wind_profile = false;
    bool terrain_aligned = false;
    if (abl) {
        pp_abl.query("initial_wind_profile", init_wind_profile);
        pp_abl.query("terrain_aligned_profile", terrain_aligned);
    }
    // ABL.initial_wind_profile only sets these fields in the interior; any
    // other scalar starts from its own initial condition, so an offset on its
    // boundary cannot conflict with the ABL profile
    const bool abl_initialized = (fld.name() == "velocity") ||
                                 (fld.name() == "temperature") ||
                                 (fld.name() == "tke");

    // Terrain lets the offset be checked against the ground it stands on,
    // using the terrain file the run reads:
    // - with TerrainDrag, its file, including its default name. TerrainDrag
    //   ignores the file and builds its terrain from the waves when OceanWaves
    //   is built before it and no vof field is (MultiPhase declares vof unless
    //   it captures the interface with a level set), so there is nothing to
    //   check against then;
    // - otherwise, the file the ABL physics reads for a terrain-aligned
    //   initial profile (TerrainDrag.terrain_file, same default).
    // Decided from the inputs rather than from the physics manager, since the
    // boundary conditions of velocity are set up before any physics is built.
    const auto terrain_drag_pos = std::ranges::find(physics, "TerrainDrag");
    const bool terrain_drag = (terrain_drag_pos != physics.end());
    const auto built_before_terrain_drag = [&](const std::string& name) {
        return std::find(physics.begin(), terrain_drag_pos, name) !=
               terrain_drag_pos;
    };
    std::string interface_model{"vof"};
    amrex::ParmParse("MultiPhase")
        .query("interface_capturing_method", interface_model);
    const bool vof_before_terrain_drag =
        built_before_terrain_drag("MultiPhase") &&
        (amrex::toLower(interface_model) != "levelset");
    const bool terrain_from_waves = terrain_drag &&
                                    built_before_terrain_drag("OceanWaves") &&
                                    !vof_before_terrain_drag;
    // Both readers default to terrain.amrwind (TerrainDrag::m_terrain_file,
    // ABLFieldInit::m_terrain_file) and stop the run when the file is missing
    std::string terrain_file;
    std::string terrain_user;
    if (terrain_drag && !terrain_from_waves) {
        terrain_user = "TerrainDrag";
    } else if (init_wind_profile && terrain_aligned) {
        terrain_user = "ABL terrain-aligned profile";
    }
    const bool terrain_required = !terrain_user.empty();
    if (terrain_required) {
        terrain_file = "terrain.amrwind";
        amrex::ParmParse("TerrainDrag").query("terrain_file", terrain_file);
    }

    // Read the terrain once for every face that is checked against it
    amrex::Vector<amrex::Real> xterrain;
    amrex::Vector<amrex::Real> yterrain;
    amrex::Vector<amrex::Real> zterrain;
    bool have_terrain = false;
    if (!terrain_file.empty()) {
        std::ifstream terrain_reader(terrain_file, std::ios::in);
        have_terrain = terrain_reader.good();
        terrain_reader.close();
        if (!have_terrain && terrain_required) {
            amrex::Abort(
                "TabulatedProfile: cannot open the " + terrain_user +
                " terrain file " + terrain_file +
                " to check the ground along the inflow faces");
        }
        if (have_terrain) {
            ioutils::read_flat_grid_file(
                terrain_file, xterrain, yterrain, zterrain);
        }
    }
    const auto& geom = fld.repo().mesh().Geom(0);
    amrex::Real ground_tol = geom.CellSize(AMREX_SPACEDIM - 1);
    pp.query("ground_tolerance", ground_tol);
    if (ground_tol < 0.0_rt) {
        amrex::Abort("TabulatedProfile.ground_tolerance must not be negative");
    }

    // The existing 1-D RANS profile file puts w in the fourth column where
    // this one puts temperature, and its readers take no header, so the two
    // cannot share a file
    std::string rans_file;
    amrex::ParmParse("ABL").query("rans_1dprofile_file", rans_file);

    const auto want = wanted_columns(fld.name(), ncomp);
    const auto& bctype = fld.bc_type();

    std::map<std::string, ProfileData> cache;
    amrex::Vector<amrex::Real> z_all;
    amrex::Vector<amrex::Real> vals_all;
    bool any_profile = false;

    for (int face = 0; face < nfaces; ++face) {
        const auto bct = bctype[face];
        if ((bct != BC::mass_inflow) && (bct != BC::mass_inflow_outflow)) {
            continue;
        }

        amrex::ParmParse pp_face(face_names[face]);
        std::string fname = default_file;
        pp_face.query("tabulated_profile_file", fname);
        // A face takes the profile only when it asks for this UDF through the
        // key of its own boundary type
        const std::string udf_key =
            fld.name() + ((bct == BC::mass_inflow) ? ".inflow_type"
                                                   : ".inflow_outflow_type");
        std::string udf_type{"ConstDirichlet"};
        pp_face.query(udf_key, udf_type);

        // One fill operator serves every inflow face of a field, so a face
        // cannot use another UDF alongside this one: it would silently get
        // the constant instead (or, the other way round, the profile would be
        // dropped). In a run the BC setup already refuses such a mix (see
        // register_velocity_dirichlet); this is kept as a defensive check for
        // direct construction.
        if ((udf_type != identifier()) && (udf_type != "ConstDirichlet")) {
            amrex::Abort(
                "TabulatedProfile: " + face_names[face] + "." + udf_key +
                " = " + udf_type +
                " cannot be combined with TabulatedProfile on another face, "
                "since one inflow UDF serves every inflow face of " +
                fld.name());
        }

        amrex::Real zoffset = default_zoffset;
        pp_face.query("tabulated_profile_zoffset", zoffset);
        m_op.zoffset[face] = zoffset;

        if ((udf_type != identifier()) || fname.empty()) {
            // Fall back to the constant value given for this face. The offset
            // and ground checks below only concern a profile, so a constant
            // face may stand on any ground.
            // Velocity defaults to zero on a face without a value, as for
            // ConstDirichlet; a scalar has no sensible default (0 K for
            // temperature), so its constant must be given
            if ((fld.name() != "velocity") && !pp_face.contains(fld.name())) {
                amrex::Abort(
                    "TabulatedProfile: " + face_names[face] +
                    " does not use the profile for " + fld.name() +
                    ", so it needs the constant value " + face_names[face] +
                    "." + fld.name());
            }
            amrex::Vector<amrex::Real> cval(ncomp, 0.0_rt);
            pp_face.queryarr(fld.name(), cval, 0, ncomp);
            for (int n = 0; n < ncomp; ++n) {
                m_op.constval[(face * AMREX_SPACEDIM) + n] = cval[n];
            }
            amrex::Print() << "TabulatedProfile: " << fld.name() << " on "
                           << face_names[face] << " uses the constant value ("
                           << ((udf_type != identifier())
                                   ? (face_names[face] + "." + udf_key +
                                      " is not TabulatedProfile")
                                   : std::string("no profile file"))
                           << ")\n";
            continue;
        }

        if (abl_initialized && init_wind_profile && !terrain_aligned &&
            (zoffset != 0.0_rt)) {
            amrex::Abort(
                "TabulatedProfile: the interior is initialized from a profile "
                "measured from the bottom of the domain, since "
                "ABL.terrain_aligned_profile is off, but the boundary on " +
                face_names[face] +
                " is offset to raised ground. Either align the initial profile "
                "with the terrain or drop the offset.");
        }
        if (have_terrain) {
            check_ground_height(
                xterrain, yterrain, zterrain, geom, face, zoffset, ground_tol);
        }

        if (!cache.contains(fname)) {
            auto prof = read_profile_file(fname);

            // The RANS profile readers take no comment lines, so a header
            // that would tell the two layouts apart would empty the RANS
            // profile; the two inputs need separate files
            if (!rans_file.empty() && same_file(fname, rans_file)) {
                amrex::Abort(
                    "TabulatedProfile: " + fname +
                    " is also used as ABL.rans_1dprofile_file, whose fourth "
                    "column is w rather than T and whose readers take no "
                    "header. Use a separate file for the inflow profile.");
            }
            if (!prof.has_header) {
                amrex::Print()
                    << "TabulatedProfile: " << fname
                    << " has no header, assuming columns z u v T"
                    << ((prof.colnames.size() > 3) ? " tke" : "") << "\n";
            }

            cache[fname] = std::move(prof);
        }
        const auto& prof = cache[fname];

        const int nz = static_cast<int>(prof.z.size());
        const int offset = static_cast<int>(z_all.size());
        m_op.offset[face] = offset;
        m_op.npts[face] = nz;
        any_profile = true;

        z_all.insert(z_all.end(), prof.z.begin(), prof.z.end());
        vals_all.resize(vals_all.size() + (static_cast<size_t>(nz) * ncomp));
        for (int n = 0; n < ncomp; ++n) {
            // A vertical velocity column is optional and defaults to zero;
            // every other component must be tabulated
            const int col = prof.column_index(want[n]);
            if ((col < 0) && (want[n] != "w")) {
                amrex::Abort(
                    "TabulatedProfile: " + fname + " has no " + want[n] +
                    " column, needed for the " + fld.name() +
                    " boundary condition on " + face_names[face] +
                    (prof.has_header
                         ? " (columns named by line " +
                               std::to_string(prof.header_line) +
                               "; if that line is a note, reword it so that "
                               "it does not start with 'z')"
                         : std::string()));
            }
            for (int k = 0; k < nz; ++k) {
                vals_all[(ncomp * (offset + k)) + n] =
                    (col < 0) ? 0.0_rt : prof.cols[col][k];
            }
        }

        if ((bct == BC::mass_inflow) && (fld.name() == "velocity")) {
            const auto& probdom = fld.repo().mesh().Geom(0).ProbDomain();
            // Heights measured from the ground of this face. On an x or y
            // face standing on raised ground the part below the ground is
            // buried in the terrain and carries no flow, so only the column
            // from the ground up is checked (none when the ground is above
            // the domain top). A bottom or top face is checked at its height.
            const bool lateral =
                (face % AMREX_SPACEDIM) != (AMREX_SPACEDIM - 1);
            const amrex::Real zbottom =
                (lateral ? amrex::max(probdom.lo(2), zoffset) : probdom.lo(2)) -
                zoffset;
            const amrex::Real ztop = probdom.hi(2) - zoffset;
            if (!lateral || (ztop > zbottom)) {
                check_inflow_direction(
                    z_all, vals_all, offset, nz, ncomp, face, zbottom, ztop,
                    fname);
            }
        }

        amrex::Print() << "TabulatedProfile: " << fld.name() << " on "
                       << face_names[face] << " from " << fname << " (" << nz
                       << " levels";
        if (zoffset != 0.0_rt) {
            amrex::Print() << ", ground at z = " << zoffset;
        }
        amrex::Print() << ")\n";
    }

    if (!any_profile) {
        amrex::Abort(
            "TabulatedProfile was requested for " + fld.name() +
            " but no inflow face uses it with a profile file. Set "
            "<face>." +
            fld.name() +
            ".inflow_type (or .inflow_outflow_type) = TabulatedProfile on an "
            "inflow face, and TabulatedProfile.filename or "
            "<face>.tabulated_profile_file.");
    }

    m_z_d.resize(z_all.size());
    m_vals_d.resize(vals_all.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, z_all.begin(), z_all.end(), m_z_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vals_all.begin(), vals_all.end(),
        m_vals_d.begin());
}

} // namespace kynema_sgf::udf
