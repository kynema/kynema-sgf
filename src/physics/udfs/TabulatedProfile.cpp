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

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

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

/** Map the spellings accepted in a header onto the field names used internally
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

/** Read the header of a profile file, if it has one
 *
 *  A leading comment line whose first entry is the height is taken as the list
 *  of column names. Any other comment line is skipped.
 */
bool parse_header(
    const std::string& line,
    amrex::Vector<std::string>& colnames,
    const std::string& fname)
{
    std::istringstream iss(line.substr(line.find('#') + 1));
    amrex::Vector<std::string> names;
    std::string name;
    while (iss >> name) {
        names.push_back(canonical_name(name));
    }
    if (names.empty() || (names[0] != "z")) {
        return false;
    }
    colnames.assign(names.begin() + 1, names.end());
    for (int i = 0; i < static_cast<int>(colnames.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(colnames.size()); ++j) {
            if (colnames[i] == colnames[j]) {
                amrex::Abort(
                    "TabulatedProfile: the header of " + fname + " names '" +
                    colnames[i] + "' more than once");
            }
        }
    }
    return true;
}

/** Where a problem was found, for a message the reader can act on
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
    return static_cast<amrex::Real>(val);
}

/** Read a whitespace-separated profile file
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
    std::string line;
    int lineno = 0;

    while (std::getline(infile, line)) {
        ++lineno;
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            continue;
        }
        if (line[first] == '#') {
            // Only a header ahead of the data names the columns
            if (rows.empty() && !prof.has_header) {
                prof.has_header = parse_header(line, prof.colnames, fname);
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
    if (prof.has_header) {
        if (static_cast<int>(prof.colnames.size()) != ncols - 1) {
            amrex::Abort(
                "TabulatedProfile: the header of " + fname + " names " +
                std::to_string(prof.colnames.size() + 1) +
                " columns but the data has " + std::to_string(ncols));
        }
    } else {
        prof.colnames = assumed_columns(ncols, fname);
    }

    const int nz = static_cast<int>(rows.size());
    prof.z.resize(nz);
    prof.cols.resize(ncols - 1);
    for (auto& col : prof.cols) {
        col.resize(nz);
    }
    for (int k = 0; k < nz; ++k) {
        prof.z[k] = rows[k][0];
        if ((k > 0) && (prof.z[k] <= prof.z[k - 1])) {
            amrex::Abort(
                "TabulatedProfile: heights must increase strictly, but " +
                at(fname, row_lines[k]) + " has " + std::to_string(prof.z[k]) +
                " after " + std::to_string(prof.z[k - 1]));
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

    // Only the part of the column the domain actually reaches matters. An
    // entry influences that range when the span between its neighbours
    // overlaps it, and the outermost entries reach beyond the table because
    // the nearest tabulated value is held outside it
    bool enters = false;
    bool leaves = false;
    for (int k = 0; k < nz; ++k) {
        const auto below =
            (k == 0) ? -constants::LARGE_NUM : heights[offset + k - 1];
        const auto above =
            (k == nz - 1) ? constants::LARGE_NUM : heights[offset + k + 1];
        if ((above < zlo) || (below > zhi)) {
            continue;
        }
        const auto un = into * vals[(ncomp * (offset + k)) + dir];
        if (un > constants::TIGHT_TOL) {
            enters = true;
        }
        if (un < -constants::TIGHT_TOL) {
            leaves = true;
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
 */
void check_ground_height(
    const std::string& terrain_file,
    const amrex::Geometry& geom,
    const int face,
    const amrex::Real zoffset,
    const amrex::Real tol)
{
    const int dir = face % AMREX_SPACEDIM;
    if (dir == AMREX_SPACEDIM - 1) {
        return;
    }

    amrex::Vector<amrex::Real> xterrain;
    amrex::Vector<amrex::Real> yterrain;
    amrex::Vector<amrex::Real> zterrain;
    ioutils::read_flat_grid_file(terrain_file, xterrain, yterrain, zterrain);

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();

    // Sample where the boundary cells sit, along the face and at its own edge
    const int tdir = (dir == 0) ? 1 : 0;
    const amrex::Real ncoord =
        (face < AMREX_SPACEDIM) ? problo[dir] : probhi[dir];
    const int npts = geom.Domain().length(tdir);

    amrex::Real zmin = constants::LARGE_NUM;
    amrex::Real zmax = -constants::LARGE_NUM;
    for (int n = 0; n < npts; ++n) {
        const auto tcoord = problo[tdir] + ((n + 0.5_rt) * dx[tdir]);
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
            " varies between " + std::to_string(zmin) + " and " +
            std::to_string(zmax) +
            ", and only a uniform lift is supported. Raise "
            "TabulatedProfile.ground_tolerance to accept this face, or drive "
            "it with a boundary plane instead.");
    }

    const auto zground = 0.5_rt * (zmin + zmax);
    if (std::abs(zground - zoffset) > tol) {
        amrex::Abort(
            "TabulatedProfile: the ground on " + face_names[face] + " is at " +
            std::to_string(zground) +
            " but the profile is offset "
            "by " +
            std::to_string(zoffset) +
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

    const int ncomp = fld.num_comp();
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

    // The interior profile and the boundary have to be measured from the same
    // place, or the inflow fights the interior at the boundary
    amrex::ParmParse pp_abl("ABL");
    bool init_wind_profile = false;
    bool terrain_aligned = false;
    pp_abl.query("initial_wind_profile", init_wind_profile);
    pp_abl.query("terrain_aligned_profile", terrain_aligned);

    // Terrain lets the offset be checked against the ground it stands on
    std::string terrain_file;
    if (amrex::ParmParse("TerrainDrag").query("terrain_file", terrain_file) ==
        0) {
        amrex::ParmParse("ImmersedTerrain").query("terrain_file", terrain_file);
    }
    const auto& geom = fld.repo().mesh().Geom(0);
    amrex::Real ground_tol = geom.CellSize(AMREX_SPACEDIM - 1);
    pp.query("ground_tolerance", ground_tol);

    // The existing 1-D RANS profile file puts w in the fourth column where
    // this one puts temperature, so the two cannot be read the same way
    std::string rans_file;
    amrex::ParmParse("ABL").query("rans_1dprofile_file", rans_file);

    std::string turbulence_model;
    amrex::ParmParse("turbulence").query("model", turbulence_model);
    const bool needs_tke = (amrex::toLower(turbulence_model) == "klaxell");

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

        amrex::Real zoffset = default_zoffset;
        pp_face.query("tabulated_profile_zoffset", zoffset);
        m_op.zoffset[face] = zoffset;

        if (init_wind_profile && !terrain_aligned && (zoffset != 0.0_rt)) {
            amrex::Abort(
                "TabulatedProfile: the interior is initialized from a profile "
                "measured from the bottom of the domain, since "
                "ABL.terrain_aligned_profile is off, but the boundary on " +
                face_names[face] +
                " is offset to raised ground. Either align the initial profile "
                "with the terrain or drop the offset.");
        }
        if (!terrain_file.empty()) {
            std::ifstream terrain_reader(terrain_file, std::ios::in);
            const bool have_terrain = terrain_reader.good();
            terrain_reader.close();
            if (have_terrain) {
                check_ground_height(
                    terrain_file, geom, face, zoffset, ground_tol);
            }
        }

        if (fname.empty()) {
            // Fall back to the constant value given for this face
            amrex::Vector<amrex::Real> cval(ncomp, 0.0_rt);
            pp_face.queryarr(fld.name(), cval, 0, ncomp);
            for (int n = 0; n < ncomp; ++n) {
                m_op.constval[(face * AMREX_SPACEDIM) + n] = cval[n];
            }
            amrex::Print() << "TabulatedProfile: " << fld.name() << " on "
                           << face_names[face]
                           << " has no profile, using the constant value\n";
            continue;
        }

        if (!cache.contains(fname)) {
            auto prof = read_profile_file(fname);

            if (!prof.has_header) {
                if (!rans_file.empty() && (fname == rans_file)) {
                    amrex::Abort(
                        "TabulatedProfile: " + fname +
                        " is also used as ABL.rans_1dprofile_file, whose "
                        "fourth column is w rather than T. Add a header line "
                        "to say which columns the file actually holds.");
                }
                amrex::Print()
                    << "TabulatedProfile: " << fname
                    << " has no header, assuming columns z u v T"
                    << ((prof.colnames.size() > 3) ? " tke" : "") << "\n";
            }

            if (needs_tke && (prof.column_index("tke") < 0)) {
                amrex::Abort(
                    "TabulatedProfile: the KLAxell model solves a TKE equation "
                    "but " +
                    fname + " has no tke column");
            }

            cache[fname] = prof;
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
                    " boundary condition on " + face_names[face]);
            }
            for (int k = 0; k < nz; ++k) {
                vals_all[(ncomp * (offset + k)) + n] =
                    (col < 0) ? 0.0_rt : prof.cols[col][k];
            }
        }

        if ((bct == BC::mass_inflow) && (fld.name() == "velocity")) {
            const auto& probdom = fld.repo().mesh().Geom(0).ProbDomain();
            check_inflow_direction(
                z_all, vals_all, offset, nz, ncomp, face,
                probdom.lo(2) - zoffset, probdom.hi(2) - zoffset, fname);
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
            " but no profile file was given. Set TabulatedProfile.filename, "
            "or <face>.tabulated_profile_file on an inflow face.");
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
