#ifndef ACTUATOR_SECTOR_OPS_H
#define ACTUATOR_SECTOR_OPS_H

#include "src/wind_energy/actuator/sector/ActuatorSector.H"
#include "src/wind_energy/actuator/actuator_ops.H"
#include "src/wind_energy/actuator/actuator_utils.H"
#include "src/core/FieldRepo.H"

#include <cmath>
#include <numbers>
#include <string>

namespace kynema_sgf::actuator::sector {

void build_radial_grid(ActuatorSectorData& meta);

vs::Tensor orientation_matrix_from_normal(const vs::Vector& normal);

vs::Tensor rotation_matrix_from_euler_degrees(const vs::Vector& angles);

vs::Tensor rotation_matrix_from_vector(const vs::Vector& rotation_vector);

vs::Tensor orientation_at_time(
    amrex::Real time,
    const vs::Tensor& initial_orientation,
    const vs::Vector& rotor_angular_velocity);

vs::Vector blade_position(
    amrex::Real radius,
    amrex::Real theta,
    amrex::Real time,
    const vs::Vector& center0,
    const vs::Vector& translation_velocity,
    const vs::Tensor& initial_orientation,
    const vs::Vector& rotor_angular_velocity);

void blade_basis(
    amrex::Real theta,
    const vs::Tensor& orientation,
    vs::Vector& e_r,
    vs::Vector& e_theta,
    vs::Vector& e_normal);

void update_midpoint_sample_points(ActuatorSector::DataType& data);

void build_gaussian_table(ActuatorSectorData& meta);

amrex::Real timestep_width(const CFDSim& sim);

void prepare_netcdf_file(
    const std::string& ncfile,
    const ActuatorSectorData& meta,
    const ActInfo& info);

void write_netcdf(
    const std::string& ncfile,
    const ActuatorSectorData& meta,
    const ActInfo& info,
    const ActGrid& grid,
    amrex::Real time,
    int time_index);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real gaussian_norm_3d()
{
    const amrex::Real pi = std::numbers::pi_v<amrex::Real>;
    return 1.0_rt / (pi * std::sqrt(pi));
}

} // namespace kynema_sgf::actuator::sector

namespace kynema_sgf::actuator::ops {

template <>
struct ReadInputsOp<ActuatorSector, ActSrcSector>
{
    void operator()(ActuatorSector::DataType& data, const utils::ActParser& pp);
};

template <>
struct InitDataOp<ActuatorSector, ActSrcSector>
{
    void operator()(ActuatorSector::DataType& data);
};

template <>
struct UpdatePosOp<ActuatorSector, ActSrcSector>
{
    void operator()(ActuatorSector::DataType& data);
};

template <>
struct UpdateVelOp<ActuatorSector, ActSrcSector>
{
    void operator()(ActuatorSector::DataType& /*data*/) {}
};

template <>
struct ComputeForceOp<ActuatorSector, ActSrcSector>
{
    void operator()(ActuatorSector::DataType& data);
};

template <>
struct ProcessOutputsOp<ActuatorSector, ActSrcSector>
{
    explicit ProcessOutputsOp(ActuatorSector::DataType& data) : m_data(data) {}
    void read_io_options(const utils::ActParser& pp);
    void prepare_outputs(const std::string& out_dir);
    void write_outputs();

private:
    ActuatorSector::DataType& m_data;
    std::string m_nc_filename;
    int m_out_freq{1};
};

template <>
class ActSrcOp<ActuatorSector, ActSrcSector>
{
private:
    ActuatorSector::DataType& m_data;
    Field& m_act_src;

    DeviceVecList m_pos;
    DeviceVecList m_force;
    DeviceVecList m_epsilon;
    amrex::Gpu::DeviceVector<amrex::Real> m_gaussian_table;

    void copy_to_device();

public:
    explicit ActSrcOp(ActuatorSector::DataType& data)
        : m_data(data)
        , m_act_src(m_data.sim().repo().get_field("actuator_src_term"))
    {}

    void initialize();

    void setup_op() { copy_to_device(); }

    void
    operator()(int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom);
};

} // namespace kynema_sgf::actuator::ops

#endif /* ACTUATOR_SECTOR_OPS_H */
