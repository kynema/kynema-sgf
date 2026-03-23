#include <cmath>

#include "amr-wind/wind_energy/ABLFieldInitFile.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

ABLFieldInitFile::ABLFieldInitFile()
{
#ifndef AMR_WIND_USE_NETCDF
    // Assert netcdf must be used for initial condition file
    amrex::Abort(
        "ABLFieldInitFile: AMR-Wind was not built with NetCDF support; "
        "initial_condition_input_file cannot be used.");
#endif
    amrex::ParmParse pp_abl("ABL");
    // Get netcdf input file name
    pp_abl.get("initial_condition_input_file", m_ic_input);
}

bool ABLFieldInitFile::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const int lev) const
{
#ifdef AMR_WIND_USE_NETCDF
    // Load the netcdf file with data if specified in the inputs
    if (lev == 0) {
        // Ensure that the input dimensions match the coarsest grid size
        const auto& domain = geom.Domain();

        // The indices that determine the start and end points of the i, j, k
        // arrays. The max and min are there to ensure that the points form
        // ghost cells are not used
        auto i0 = amrex::max(vbx.smallEnd(0), domain.smallEnd(0));
        auto i1 = amrex::min(vbx.bigEnd(0), domain.bigEnd(0));

        auto j0 = amrex::max(vbx.smallEnd(1), domain.smallEnd(1));
        auto j1 = amrex::min(vbx.bigEnd(1), domain.bigEnd(1));

        auto k0 = amrex::max(vbx.smallEnd(2), domain.smallEnd(2));
        auto k1 = amrex::min(vbx.bigEnd(2), domain.bigEnd(2));

        // Loop through all points in the domain and set velocities to values
        // from the input file
        // start is the first index from where to read data
        amrex::Vector<size_t> start{
            static_cast<size_t>(i0), static_cast<size_t>(j0),
            static_cast<size_t>(k0)};
        // count is the total number of elements to read in each direction
        amrex::Vector<size_t> count{
            static_cast<size_t>(i1 - i0 + 1), static_cast<size_t>(j1 - j0 + 1),
            static_cast<size_t>(k1 - k0 + 1)};

        // Working vector to read data onto host
        const auto dlen = count[0] * count[1] * count[2];
        // Vector to store the 3d data into a single array and set size
        amrex::Gpu::DeviceVector<amrex::Real> uvel_d(dlen, 0.0_rt);
        amrex::Gpu::DeviceVector<amrex::Real> vvel_d(dlen, 0.0_rt);
        amrex::Gpu::DeviceVector<amrex::Real> wvel_d(dlen, 0.0_rt);

        // Open the netcdf input file
        // This file should have the same dimensions as the simulation
#ifdef AMREX_USE_OMP
#pragma omp critical(ncf_open_close)
#endif
        {
            auto ncf = ncutils::NCFile::open(m_ic_input, NC_NOWRITE);

            // The x, y and z velocity components (u, v, w)
            auto uvel = ncf.var("uvel");
            auto vvel = ncf.var("vvel");
            auto wvel = ncf.var("wvel");

            // Read the velocity components u, v, w and copy to device
            amrex::Vector<amrex::Real> tmp(dlen, 0.0_rt);
            uvel.get(tmp.data(), start, count);
            amrex::Gpu::copy(
                amrex::Gpu::hostToDevice, tmp.begin(), tmp.end(),
                uvel_d.begin());
            vvel.get(tmp.data(), start, count);
            amrex::Gpu::copy(
                amrex::Gpu::hostToDevice, tmp.begin(), tmp.end(),
                vvel_d.begin());
            wvel.get(tmp.data(), start, count);
            amrex::Gpu::copy(
                amrex::Gpu::hostToDevice, tmp.begin(), tmp.end(),
                wvel_d.begin());

            // Close the netcdf file
            ncf.close();
        }

        // Pointers to velocity objects
        const auto* uvel_dptr = uvel_d.data();
        const auto* vvel_dptr = vvel_d.data();
        const auto* wvel_dptr = wvel_d.data();

        // Get count components for device
        const auto ct1 = count[1];
        const auto ct2 = count[2];

        // Amrex parallel for to assign the velocity at each point
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // The counter to go from 3d to 1d vector
            auto idx = (i - i0) * ct2 * ct1 + (j - j0) * ct2 + (k - k0);
            // Pass values from temporary array to the velocity field
            velocity(i, j, k, 0) = uvel_dptr[idx];
            velocity(i, j, k, 1) = vvel_dptr[idx];
            velocity(i, j, k, 2) = wvel_dptr[idx];
        });
        // Populated directly, do not fill from another level
        return false;
    } // Skip level and interpolate data from already loaded coarse levels
    return true;

#else
    amrex::ignore_unused(vbx, geom, velocity, lev);
    return false;
#endif
}

} // namespace amr_wind
