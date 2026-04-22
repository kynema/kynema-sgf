#ifndef UNIFORM_CT_OPS_H
#define UNIFORM_CT_OPS_H

#include <numbers>
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "amr-wind/wind_energy/actuator/disk/UniformCt.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
namespace amr_wind::actuator::ops {
namespace uniformct {

void parse_and_gather_params(const utils::ActParser& pp, UniformCtData& data);
void update_disk_points(UniformCt::DataType& data);
} // namespace uniformct

template <>
struct ReadInputsOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data, const utils::ActParser& pp)
    {
        auto& meta = data.meta();
        uniformct::parse_and_gather_params(pp, meta);
        base::do_parse_based_computations<UniformCt>(data);
    }
};

template <>
struct UpdatePosOp<UniformCt, ActSrcDisk>
{
    void operator()(typename UniformCt::DataType& data)
    {
        if (!data.sim().helics().is_activated()) {
            return;
        }

#ifdef AMR_WIND_USE_HELICS

        const amrex::Real wind_direction =
            -data.sim().helics().m_turbine_yaw_to_amrwind[data.info().id] +
            90.0_rt;
        const amrex::Real wind_direction_radian =
            amr_wind::utils::radians(wind_direction);

        auto& meta = data.meta();

        meta.normal_vec[0] = std::cos(wind_direction_radian);
        meta.normal_vec[1] = std::sin(wind_direction_radian);
        meta.normal_vec[2] = 0.0_rt;

        meta.sample_vec[0] = meta.normal_vec[0];
        meta.sample_vec[1] = meta.normal_vec[1];
        meta.sample_vec[2] = meta.normal_vec[2];

        if (data.info().is_root_proc) {
            std::cout << "turbine " << data.info().id
                      << " normal: " << meta.normal_vec[0] << ' '
                      << meta.normal_vec[1] << std::endl;
        }

        uniformct::update_disk_points(data);
#endif
    }
};
// compute the total force over a given radial section
// this will then be spread uniformly over that area in the
// projection step
template <>
struct ComputeForceOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data)
    {
        auto& grid = data.grid();
        auto& ddata = data.meta();
        const auto& normal = ddata.normal_vec;
        const amrex::Real uInfSqr = base::compute_reference_velocity_sqr(ddata);
        base::set_thrust_coefficient(ddata, uInfSqr);

        const int npts = ddata.num_force_pts;
        const amrex::Real rho = ddata.density;
        const amrex::Real aeroPressure =
            0.5_rt * uInfSqr * rho * ddata.current_ct;
        const amrex::Real dr = ddata.dr;

        for (int ip = 0; ip < npts; ++ip) {
            const amrex::Real r = (ip + 0.5_rt) * dr;
            const amrex::Real rp = r + (dr * 0.5_rt);
            const amrex::Real rm = r - (dr * 0.5_rt);
            const amrex::Real a =
                std::numbers::pi_v<amrex::Real> * (rp * rp - rm * rm);
            // disk force should always point in direction of the normal
            grid.force[ip] = (aeroPressure * a) * normal;
        }
#ifdef AMR_WIND_USE_HELICS

        if (data.info().is_root_proc && data.sim().helics().is_activated()) {

            const amrex::Real ct = ddata.current_ct;
            std::cout << "turbine" << data.info().id << " thrust: " << ct
                      << std::endl;
            const amrex::Real cp = 0.5_rt * (ct + ct * std::sqrt(1.0_rt - ct));
            const amrex::Real area = std::numbers::pi_v<amrex::Real> * 0.25_rt *
                                     ddata.diameter * ddata.diameter;
            const amrex::Real power =
                cp * 0.5_rt * rho * std::pow(uInfSqr, 1.5_rt) * area;

            data.sim().helics().m_turbine_power_to_controller[data.info().id] =
                power;
            const amrex::Real turbine_angle = std::atan2(
                ddata.reference_velocity[1], ddata.reference_velocity[0]);
            data.sim()
                .helics()
                .m_turbine_wind_direction_to_controller[data.info().id] =
                -amr_wind::utils::degrees(turbine_angle) + 270.0_rt;
        }
#endif
    }
};

template <>
struct InitDataOp<UniformCt, ActSrcDisk>
{
    void operator()(typename UniformCt::DataType& data)
    {
        ops::base::allocate_basic_grid_quantities<UniformCt>(data);
        uniformct::update_disk_points(data);
    }
};
template <>
struct ProcessOutputsOp<UniformCt, ActSrcDisk>
{
private:
    // cppcheck-suppress uninitMemberVarPrivate
    const UniformCt::DataType& m_data;
    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir;

    //! NetCDF output filename for this turbine
    std::string m_nc_filename;

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp<UniformCt, ActSrcDisk>(
        const UniformCt::DataType& data)
        : m_data(data)
    {}
    void operator()(UniformCt::DataType& /*unused*/) {}
    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }
    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        disk::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }
    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        disk::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};
} // namespace amr_wind::actuator::ops

#endif /* UNIFORM_CT_OPS_H */
