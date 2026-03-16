#ifndef JOUKOWSKY_H
#define JOUKOWSKY_H

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

struct JoukowskyData : public DiskBaseData
{
    RealList angular_velocity;
    RealList tip_correction;
    RealList root_correction;
    RealList f_normal;
    RealList f_theta;
    amrex::Real current_angular_velocity{0.0_rt};
    amrex::Real current_tip_speed_ratio{0.0_rt};
    amrex::Real vortex_core_size;
    amrex::Real current_cp;
    amrex::Real current_power;
    // --- Sorensen 2020 equation 10 constants ----
    amrex::Real root_correction_coefficient{2.335_rt};
    amrex::Real root_correction_exponent{4.0_rt};
    // --- Sorensen 2022 equation 26 constants ----
    // Note that Ct_rated is the Ct in the region 2 part of the power curve.
    // Ct_rated is called ct_region2 in the input file.
    amrex::Real Ct_rated{0.0_rt};
    // Note that the defaults for S0 have changed from the paper.
    // Sorensen 2022 had S0_alpha1 = 0.08_rt and S0_alpha2 = 0.05_rt.
    // Current version uses negative alpha coefficients, possibly to account
    // for the change in rotor rotation.
    amrex::Real S0_alpha1{-0.20_rt};
    amrex::Real S0_alpha2{-0.05_rt};
    // --------------------------------------------
    int num_blades{3};
    vs::Vector disk_force{0.0_rt, 0.0_rt, 0.0_rt};
    bool use_tip_correction{true};
    bool use_root_correction{true};
    amrex::Real awc_amplitude{0.0_rt};
    amrex::Real awc_angular_frequency{0.0_rt};
    int awc_azimuthal_mode{0};
    amrex::Real awc_clocking_angle{0.0_rt};
};

struct Joukowsky : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = JoukowskyData;
    using DataType = ActDataHolder<Joukowsky>;

    static std::string identifier() { return "Joukowsky"; }
};

} // namespace amr_wind::actuator

#endif /* Joukowsky_H */
