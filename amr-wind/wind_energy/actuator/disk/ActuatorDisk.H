#ifndef ACTUATORDISK_H
#define ACTUATORDISK_H

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

struct DiskBaseData
{
    int num_force_pts;
    int num_force_theta_pts{1};
    int num_vel_pts{1};
    int num_vel_pts_r{1};
    int num_vel_pts_t{1};
    vs::Vector center;
    vs::Vector normal_vec{0, 1, 0};
    vs::Vector coplanar_vec{1, 0, 0};
    vs::Vector sample_vec{0, 1, 0};
    vs::Vector reference_velocity;
    vs::Vector mean_disk_velocity;
    amrex::Real epsilon;
    amrex::Real diameter;
    amrex::Real dr;
    amrex::Real density;
    amrex::Real current_ct;
    amrex::Real diameters_to_sample{0.0_rt};
    RealList thrust_coeff;
    RealList table_velocity;
    std::string spreading_type{"LinearBasis"};
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real radius() const
    {
        return 0.5_rt * diameter;
    }
};

} // namespace amr_wind::actuator

#endif /* ACTUATORDISK_H */
