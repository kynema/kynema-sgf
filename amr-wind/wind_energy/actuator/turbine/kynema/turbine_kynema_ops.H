#ifndef TURBINE_KYNEMA_OPS_H
#define TURBINE_KYNEMA_OPS_H

#include "amr-wind/wind_energy/actuator/turbine/kynema/TurbineKynema.H"
#include "amr-wind/wind_energy/actuator/turbine/external/turbine_external_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops {

template <typename SrcTrait>
struct ReadInputsOp<TurbineKynema, SrcTrait>
{
    void operator()(TurbineKynema::DataType& data, const utils::ActParser& pp)
    {
        auto& tf = data.meta().ext_data;
        // Only one mode available for now
        const std::string sim_mode{"init"};

        if (sim_mode == "init") {
            tf.sim_mode = ::ext_turb::SimMode::init;
            amrex::Print() << "Initializing turbine:" << tf.tlabel << std::endl;
        } else if (sim_mode == "replay") {
            tf.sim_mode = ::ext_turb::SimMode::replay;
            amrex::Print() << "Replaying turbine:" << tf.tlabel << std::endl;
        } else if (sim_mode == "restart") {
            tf.sim_mode = ::ext_turb::SimMode::restart;
            amrex::Print() << "Restarting turbine:" << tf.tlabel << std::endl;
        } else {
            amrex::Abort(
                "Actuator: Invalid Kynema simulation mode: " + sim_mode);
        }

        if (tf.sim_mode == ::ext_turb::SimMode::restart) {
            pp.get("kynema_restart_file", tf.checkpoint_file);
        } else {
            pp.get("kynema_input_file", tf.input_file);
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            amrex::FileSystem::Exists(tf.input_file),
            "turbine_kynema_ops: Cannot find Kynema input file: " +
                tf.input_file);

        pp.query(
            "controller_shared_library_path", tf.controller_shared_lib_path);
        pp.query("controller_input_file", tf.controller_input_file);

        // Get hub height and turbine diameter from turbine input,
        // not AMR-Wind input file
        const YAML::Node wio = YAML::LoadFile(tf.input_file);
        {
            const auto& wio_blade = wio["components"]["blade"];
            const auto& ref_axis = wio_blade["reference_axis"];
            const auto z_values =
                ref_axis["z"]["values"].as<std::vector<amrex::Real>>();
            data.info().rotor_diameter = 2.0_rt * z_values[z_values.size() - 1];
        }
        {
            const auto& wio_tower = wio["components"]["tower"];
            const auto& tower_ref_axis = wio_tower["reference_axis"];
            const auto z_values =
                tower_ref_axis["z"]["values"].as<std::vector<amrex::Real>>();
            data.info().hub_height = z_values[z_values.size() - 1];
        }

        read_ops<TurbineKynema::DataType>(data, pp);

        // Do AMR-Wind inputs specific to kynema turbine
        tf.fluid_density = data.meta().density;

        tf.dt_ext = -1.;
        pp.query("dt", tf.dt_ext);
        if (tf.dt_ext < 0.0_rt) {
            tf.dt_ext = tf.dt_cfd;
        }

        // Number of beam elements read in from file
        pp.get("num_struct_nodes_blade", tf.num_blade_elem);
        // Initial or constant angular rotation speed
        pp.query("rot_speed_radps", tf.rotational_speed);
        if (!pp.contains("rot_speed_radps")) {
            pp.query("rot_speed_rpm", tf.rotational_speed);
            tf.rotational_speed *=
                2.0_rt * std::numbers::pi_v<amrex::Real> / 60.0_rt;
        }
        // Initial generator speed
        pp.query("generator_power_init", tf.generator_power);
        // Initial wind speed
        amrex::Vector<amrex::Real> vel_vec{tf.wind_speed, 0., 0.};
        pp.queryarr("hub_wind_vector_init", vel_vec);
        tf.wind_speed = std::sqrt(
            vel_vec[0] * vel_vec[0] + vel_vec[1] * vel_vec[1] +
            vel_vec[2] * vel_vec[2]);
        // Initial or constant nacelle yaw
        pp.query("yaw_rad", tf.yaw);
        if (!pp.contains("yaw_rad")) {
            pp.query("yaw_deg", tf.yaw);
            tf.yaw *= std::numbers::pi_v<amrex::Real> / 180.0_rt;
        }
        // Generator efficiency
        pp.query("generator_efficiency", tf.generator_efficiency);
    }
};

template <>
void determine_influenced_procs<TurbineKynema>(
    typename TurbineKynema::DataType& data)
{
    external_determine_influenced_procs<TurbineKynema::DataType>(data);
}

template <>
void determine_root_proc<TurbineKynema>(
    typename TurbineKynema::DataType& data, amrex::Vector<int>& act_proc_count)
{
    external_determine_root_proc<
        TurbineKynema::DataType, ::ext_turb::KynemaTurbine,
        ::ext_turb::KynemaSolverData>(data, act_proc_count);
}

template <typename SrcTrait>
struct InitDataOp<TurbineKynema, SrcTrait>
{
    void operator()(TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::InitDataOp<TurbineKynema>");
        init_data_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct UpdatePosOp<TurbineKynema, SrcTrait>
{
    void operator()(typename TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::UpdatePosOp<TurbineKynema>");
        update_pos_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct UpdateVelOp<TurbineKynema, SrcTrait>
{
    void operator()(typename TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::UpdateVelOp<TurbineKynema>");
        update_vel_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct ComputeForceOp<TurbineKynema, SrcTrait>
{
    void operator()(typename TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::ComputeForceOp<TurbineKynema>");
        compute_force_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct ProcessOutputsOp<TurbineKynema, SrcTrait>
{
private:
    typename TurbineKynema::DataType& m_data;

    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir;

    //! NetCDF output filename for this turbine
    std::string m_nc_filename;

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename TurbineKynema::DataType& data)
        : m_data(data)
    {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        utils::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if ((m_out_freq > 0) && (tidx % m_out_freq != 0)) {
            return;
        }

        utils::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace amr_wind::actuator::ops

#endif /* TURBINE_KYNEMA_OPS_H */
