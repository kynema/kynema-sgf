#ifndef DRONE_OPS_H
#define DRONE_OPS_H

#include "src/wind_energy/actuator/drone/Drone.H"
#include "src/wind_energy/actuator/actuator_ops.H"

namespace kynema_sgf::actuator::ops {

template <>
struct ReadInputsOp<Drone, ActSrcDrone>
{
    void operator()(Drone::DataType& data, const utils::ActParser& pp);
};

template <>
struct InitDataOp<Drone, ActSrcDrone>
{
    void operator()(Drone::DataType& data);
};

template <>
struct UpdatePosOp<Drone, ActSrcDrone>
{
    void operator()(Drone::DataType& data);
};

template <>
struct UpdateVelOp<Drone, ActSrcDrone>
{
    void operator()(Drone::DataType& data);
};

template <>
struct ComputeForceOp<Drone, ActSrcDrone>
{
    void operator()(Drone::DataType& data);
};

template <>
struct RefinementGeometryOp<Drone, ActSrcDrone>
{
    [[nodiscard]] amrex::Vector<ActuatorRefinementGeometry>
    operator()(const Drone::DataType& data, amrex::Real time) const;
};

template <>
struct ProcessOutputsOp<Drone, ActSrcDrone>
{
    explicit ProcessOutputsOp(Drone::DataType& data) : m_data(data) {}
    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_output_frequency);
    }
    void prepare_outputs(const std::string& out_dir);
    void write_outputs();

private:
    Drone::DataType& m_data;
    std::string m_nc_filename;
    int m_output_frequency{1};
};

template <>
class ActSrcOp<Drone, ActSrcDrone>
{
public:
    explicit ActSrcOp(Drone::DataType& data) : m_data(data) {}
    void initialize();
    void setup_op();
    void
    operator()(int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom);

private:
    Drone::DataType& m_data;
};

} // namespace kynema_sgf::actuator::ops

#endif
