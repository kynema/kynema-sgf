#ifndef DISK_OPS_H
#define DISK_OPS_H

#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {
namespace disk {
// refactor these when we add more disk types later
void prepare_netcdf_file(
    const std::string& name,
    const DiskBaseData& data,
    const ActInfo& info,
    const ActGrid& grid);

void write_netcdf(
    const std::string& name,
    const DiskBaseData& data,
    const ActInfo& info,
    const ActGrid& /*unused*/,
    amrex::Real time);
} // namespace disk

namespace ops {
namespace base {

class AreaComputer
{
public:
    AreaComputer(amrex::Real radius, int num_r, int num_theta);

    [[nodiscard]] amrex::Real area_section(int iRadius) const;

    [[nodiscard]] amrex::Real weight(int iRadius) const;

private:
    const amrex::Real m_area;
    const amrex::Real m_geometry_factor;
};

vs::Vector compute_coplanar_vector(const vs::Vector& normal);

void collect_parse_conflicts(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss);
void collect_parse_dependencies_one_way(
    const utils::ActParser& pp,
    const std::string& independent,
    const std::string& dependent,
    std::ostringstream& ss);
void collect_parse_dependencies(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss);

void check_error_stream(const std::ostringstream& error_collector);

void required_parameters(DiskBaseData& meta, const utils::ActParser& pp);

void optional_parameters(DiskBaseData& meta, const utils::ActParser& pp);

std::ostringstream check_for_parse_conflicts(const utils::ActParser& pp);

void compute_and_normalize_coplanar_vector(DiskBaseData& meta);

void final_checks(const DiskBaseData& meta);

amrex::RealBox compute_bounding_box(const DiskBaseData& meta);

template <typename T>
void allocate_basic_grid_quantities(typename T::DataType& data)
{
    auto& grid = data.grid();
    auto& meta = data.meta();
    // only resize the members we are going to use
    grid.pos.resize(meta.num_force_pts);
    grid.force.resize(meta.num_force_pts);
    grid.vel.resize(meta.num_vel_pts);
    grid.vel_pos.resize(meta.num_vel_pts);
    grid.density.resize(meta.num_vel_pts);
}

template <typename T>
void do_parse_based_computations(ActDataHolder<T>& data)
{
    auto& meta = data.meta();
    auto& info = data.info();

    compute_and_normalize_coplanar_vector(meta);
    info.bound_box = compute_bounding_box(meta);
}

void compute_disk_points(
    const DiskBaseData& meta,
    VecList& points,
    const vs::Vector& cylAxis,
    int offset,
    amrex::Real dOffset);

template <typename T>
void parse_and_gather_params(const utils::ActParser& pp, T& metaData)
{
    check_for_parse_conflicts(pp);
    required_parameters(metaData, pp);
    optional_parameters(metaData, pp);
}

AMREX_INLINE amrex::Real compute_reference_velocity_sqr(DiskBaseData& data)
{
    const auto& normal = data.normal_vec;
    // only use velocity components normal to disk
    // make sure uInfSqr is always positive by squaring the projection of
    // reference_velocity onto the normal
    return (data.reference_velocity & normal) *
           (data.reference_velocity & normal);
}

// we use a setter function since we want to cache this variable for netcdf
// output
AMREX_INLINE void
set_thrust_coefficient(DiskBaseData& data, const amrex::Real& uInfSqr)
{
    const amrex::Real uInfMag = std::sqrt(uInfSqr);
    data.current_ct = ::amr_wind::interp::linear(
        data.table_velocity, data.thrust_coeff, uInfMag);
}

} // namespace base

template <typename ActTrait>
struct UpdateVelOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<DiskType, ActTrait>>>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& meta = data.meta();
        const auto& grid = data.grid();
        meta.reference_velocity = {0.0_rt, 0.0_rt, 0.0_rt};
        meta.mean_disk_velocity = {0.0_rt, 0.0_rt, 0.0_rt};
        meta.density = 0.0_rt;
        auto& refVel = meta.reference_velocity;
        auto& diskVel = meta.mean_disk_velocity;
        auto& density = meta.density;
        const int np = meta.num_vel_pts / 2;
        const int num_r = meta.num_vel_pts_r;
        const int num_t = meta.num_vel_pts_t;

        const base::AreaComputer area(meta.radius(), num_r, num_t);

        for (int i = 0; i < num_r; i++) {
            const auto weight = area.weight(i);
            for (int j = 0; j < num_t; ++j) {
                density += grid.density[j + (i * num_t)] * weight;
                refVel = refVel + (grid.vel[j + (i * num_t)] * weight);
                diskVel = diskVel + (grid.vel[j + (i * num_t) + np] * weight);
            }
        }
    }
};

template <typename ActTrait>
struct UpdatePosOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<DiskType, ActTrait>>>
{
    void operator()(typename ActTrait::DataType& /*data*/) {}
};
} // namespace ops
} // namespace amr_wind::actuator

#endif /* DISK_OPS_H */
