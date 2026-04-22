#include <fstream>

#include "amr-wind/utilities/io_utils.H"

namespace amr_wind::ioutils {

void goto_next_line(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void read_flat_grid_file(
    const std::string& fname,
    amrex::Vector<amrex::Real>& xs,
    amrex::Vector<amrex::Real>& ys,
    amrex::Vector<amrex::Real>& zs)
{
    std::ifstream file(fname, std::ios::in);

    if (!file.good()) {
        amrex::Abort("Cannot find file");
    }

    size_t nx = 0;
    size_t ny = 0;
    if (!(file >> nx)) {
        amrex::Abort("Failed to read grid dimension nx");
    }
    if (!(file >> ny)) {
        amrex::Abort("Failed to read grid dimension ny");
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nx > 0, "nx must be > 0");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ny > 0, "ny must be > 0");

    xs.resize(nx);
    ys.resize(ny);
    zs.resize(nx * ny);

    for (size_t n = 0; n < nx; n++) {
        if (!(file >> xs[n])) {
            amrex::Abort("Failed to read xs[" + std::to_string(n) + "]");
        }
    }
    for (size_t n = 0; n < ny; n++) {
        if (!(file >> ys[n])) {
            amrex::Abort("Failed to read ys[" + std::to_string(n) + "]");
        }
    }
    for (size_t n = 0; n < nx * ny; n++) {
        if (!(file >> zs[n])) {
            amrex::Abort("Failed to read zs[" + std::to_string(n) + "]");
        }
    }

    file.close();
}
} // namespace amr_wind::ioutils
