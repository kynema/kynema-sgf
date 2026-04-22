#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <chrono>
#include <ctime>
#include <string>
#include <sstream>
#include <unordered_set>
#include "AMReX_Vector.H"

namespace amrex {
const char* buildInfoGetGitHash(int i);
}

namespace amr_wind::ioutils {

void goto_next_line(std::istream& is);

AMREX_INLINE std::string amr_wind_version()
{
    return {amrex::buildInfoGetGitHash(1)};
}

AMREX_INLINE std::string timestamp()
{
    auto etime = std::chrono::system_clock::now();
    auto etimet = std::chrono::system_clock::to_time_t(etime);
    amrex::Array<char, 64> buf;
    ctime_r(&etimet, buf.begin());
    std::string tstamp(buf.begin());
    return tstamp;
}

/** Generate field names for multi-component fields
 *
 *  Function generates variable names based on the number of components in the
 *  field. If the field is a scalar, it just inserts the variable name. For
 *  vectors it appends "x, y, z" to the components, and for all other fields it
 *  numbers the components start with 0.
 */
AMREX_INLINE void add_var_names(
    amrex::Vector<std::string>& vnames,
    const std::string& fname,
    const int ncomp)
{
    const amrex::Vector<std::string> comp{"x", "y", "z"};

    switch (ncomp) {
    case 1:
        vnames.push_back(fname);
        break;

    case AMREX_SPACEDIM:
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            vnames.push_back(fname + comp[i]);
        }
        break;

    default:
        for (int i = 0; i < ncomp; ++i) {
            vnames.push_back(fname + std::to_string(i));
        }
    }
}

template <typename T>
AMREX_INLINE bool all_distinct(const amrex::Vector<T>& vec)
{
    std::unordered_set<T> s(vec.begin(), vec.end());
    return static_cast<amrex::Long>(s.size()) == vec.size();
}

AMREX_INLINE void assert_with_message(const bool val, const std::string& msg)
{
    if (!val) {
        amrex::Abort(msg);
    }
}

/** Read a flattened 2D grid file
 *
 *  File is assumed to contain a single column with the following data:
 *  [nx, ny, x values (len = nx), y values (len = ny), z values (len = nx * ny)]
 */
void read_flat_grid_file(
    const std::string& fname,
    amrex::Vector<amrex::Real>& xs,
    amrex::Vector<amrex::Real>& ys,
    amrex::Vector<amrex::Real>& zs);

} // namespace amr_wind::ioutils

#endif /* IO_UTILS_H */
