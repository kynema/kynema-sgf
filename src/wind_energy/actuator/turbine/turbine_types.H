#ifndef TURBINE_TYPES_H
#define TURBINE_TYPES_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

struct FLLCData;

struct TurbineInfo : public ActInfo
{
    //! Base location of this turbine
    vs::Vector base_pos;

    //! Rotor diameter for this turbine
    amrex::Real rotor_diameter{0.0_rt};

    //! Hub height of this turbine
    amrex::Real hub_height{0.0_rt};

    TurbineInfo(const std::string& label_in, const int id_in)
        : ActInfo(label_in, id_in)
    {}
};

using DeviceVecComponent = amrex::Gpu::DeviceVector<ComponentView>;
struct TurbineBaseData
{
    //! Identifier string for turbine solver type
    std::string solver_name{"Internal"};

    //! Number of blades
    int num_blades{3};

    //! Number of actuator nodes per blade
    int num_pts_blade;

    //! Number of velocity nodes per blade
    int num_vel_pts_blade;

    //! Number of actuator nodes for tower
    int num_pts_tower;

    RealList chord;
    VecList vel_rel;

    //! Center of rotation for the rotor
    vs::Vector rot_center;

    /** Reference frame for the rotor plane
     *
     *  x - Pointing downwind along shaft/rotation axis
     *  y - lateral direction
     *  z - normal direction (pointing mostly along z-direction)
     */
    vs::Tensor rotor_frame;

    //! General epsilon provided
    vs::Vector eps_inp;

    //! Gaussian smearing factor for the blade
    vs::Vector eps_chord;

    //! Minimum epsilon when using chord based epsilon
    vs::Vector eps_min;

    //! Gaussian smearing factor for tower
    vs::Vector eps_tower;

    //! Drag coefficient for nacelle
    amrex::Real nacelle_cd{0.0_rt};

    //! Wetted surface area for nacelle
    amrex::Real nacelle_area{0.0_rt};

    std::vector<ComponentView> blades;
    std::vector<FLLCData> fllc;
    ComponentView tower;
    ComponentView hub;
};

struct TurbineType : public ActuatorType
{};

} // namespace amr_wind::actuator

#endif /* TURBINE_TYPES_H */
