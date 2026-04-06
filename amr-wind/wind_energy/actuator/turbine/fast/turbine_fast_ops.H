#ifndef TURBINE_FAST_OPS_H
#define TURBINE_FAST_OPS_H

#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.H"
#include "amr-wind/wind_energy/actuator/turbine/external/turbine_external_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops {

template <typename SrcTrait>
struct ReadInputsOp<TurbineFast, SrcTrait>
{
    void operator()(TurbineFast::DataType& data, const utils::ActParser& pp)
    {
        read_ops<TurbineFast::DataType>(data, pp);

        auto& tf = data.meta().ext_data;

        pp.get("openfast_start_time", tf.start_time);
        pp.get("openfast_stop_time", tf.stop_time);

        std::string sim_mode = (tf.start_time > 0.0_rt) ? "restart" : "init";
        pp.query("openfast_sim_mode", sim_mode);

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
                "Actuator: Invalid OpenFAST simulation mode: " + sim_mode);
        }

        // If we are using OpenFAST restart file, require that the user provide
        // the path to the checkpoint file.
        if (tf.sim_mode == ::ext_turb::SimMode::restart) {
            pp.get("openfast_restart_file", tf.checkpoint_file);
        } else {
            pp.get("openfast_input_file", tf.input_file);
        }

        // Check if user wants to override density checks
        bool override_density_check = false;
        pp.query("override_density_check", override_density_check);
        // Check if simulation is a restart
        bool is_restart = data.sim().io_manager().is_restart();

        perform_checks(data);

        perform_density_checks(
            data.meta().density, override_density_check, is_restart);
    }

    void perform_checks(typename TurbineFast::DataType& data)
    {
        const auto& time = data.sim().time();
        // Ensure that we are using fixed timestepping scheme
        AMREX_ALWAYS_ASSERT(!time.adaptive_timestep());
    }

    // Check other density parsing arguments against turbine air density
    void perform_density_checks(
        const amrex::Real rho_tdef, const bool no_check, const bool is_rst)
    {
        // Check if other arguments are in input file
        amrex::ParmParse pp_incflo("incflo");
        amrex::ParmParse pp_cstdns("ConstValue.density");
        bool is_incflo = pp_incflo.contains("density");
        bool is_cstdns = pp_cstdns.contains("value");
        // Flags for detecting conflicts
        bool cft_incflo = false;
        bool cft_cstdns = false;

        // Do warnings or aborts depending on values
        const amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();
        if (is_incflo) {
            amrex::Real rho_incflo = 1.0_rt;
            pp_incflo.query("density", rho_incflo);
            if (std::abs(rho_incflo - rho_tdef) > tiny) {
                cft_incflo = true;
                if (!no_check) {
                    amrex::Abort(
                        "Density conflict detected between TurbineFast "
                        "density and incflo.density.\n----- Values are " +
                        std::to_string(rho_tdef) + " and " +
                        std::to_string(rho_incflo) +
                        ", respectively. Check the problem setup and the "
                        "turbine definition.\n----- If this difference is "
                        "intended, set the override_density_check flag to "
                        "true.\n");
                } else {
                    amrex::Print()
                        << "WARNING: Density conflict detected between FAST "
                           "Turbine density and incflo.density.\n----- Values "
                           "are " +
                               std::to_string(rho_tdef) + " and " +
                               std::to_string(rho_incflo) +
                               ", respectively. Abort overridden by flag.\n";
                }
            }
        }
        if (is_cstdns) {
            amrex::Real rho_cstdns = 1.0_rt;
            pp_cstdns.query("value", rho_cstdns);
            if (std::abs(rho_cstdns - rho_tdef) > tiny) {
                cft_cstdns = true;
                if (!no_check) {
                    amrex::Abort(
                        "Density conflict detected between TurbineFast "
                        "density and ConstValue.density.value.\n----- Values "
                        "are " +
                        std::to_string(rho_tdef) + " and " +
                        std::to_string(rho_cstdns) +
                        ", respectively. Check the problem setup and the "
                        "turbine definition.\n----- If this difference is "
                        "intended, set the override_density_check flag to "
                        "true.\n");
                } else {
                    amrex::Print()
                        << "WARNING: Density conflict detected between "
                           "TurbineFast "
                           "density and "
                           "ConstValue.density.value.\n----- Values are " +
                               std::to_string(rho_tdef) + " and " +
                               std::to_string(rho_cstdns) +
                               ", respectively. Abort overridden by flag.\n";
                }
            }
        }
        if (!is_incflo && !is_cstdns) {
            amrex::Print()
                << "TurbineFast: no density conflict detected, but no density "
                   "arguments found to check against.\n";
        } else if (!cft_incflo && !cft_cstdns) {
            amrex::Print() << "TurbineFast: no density conflict detected.\n";
        }
        if (is_rst) {
            amrex::Print()
                << "TurbineFast: simulation begins using a checkpoint file, "
                   "density compatibility must be manually confirmed between "
                   "the precursor simulation and the specified TurbineFast "
                   "density, " +
                       std::to_string(rho_tdef)
                << ".\n";
        }
    }
};

template <>
inline void
determine_influenced_procs<TurbineFast>(typename TurbineFast::DataType& data)
{
    external_determine_influenced_procs<TurbineFast::DataType>(data);
}

template <>
inline void determine_root_proc<TurbineFast>(
    typename TurbineFast::DataType& data, amrex::Vector<int>& act_proc_count)
{
    external_determine_root_proc<
        TurbineFast::DataType, ::ext_turb::FastTurbine,
        ::ext_turb::FastSolverData>(data, act_proc_count);
}

template <typename SrcTrait>
struct InitDataOp<TurbineFast, SrcTrait>
{
    void operator()(TurbineFast::DataType& data)
    {
        BL_PROFILE("amr-wind::InitDataOp<TurbineFast>");
        // Ensure that FAST simulation time is set properly before doing any
        // initialization tasks. We perform this check here to account for
        // restart which is only known after reading the checkpoint file.
        check_fast_sim_time(data);

        init_data_op<TurbineFast::DataType>(data);
    }

    void check_fast_sim_time(TurbineFast::DataType& data)
    {
        const auto& time = data.sim().time();

        // Set OpenFAST solver end time to be at least as long as the CFD time.
        // User can choose a longer duration in input file.
        const amrex::Real stop1 = time.stop_time() > 0.0_rt
                                      ? time.stop_time()
                                      : std::numeric_limits<amrex::Real>::max();
        const amrex::Real stop2 = time.stop_time_index() > -1
                                      ? time.stop_time_index() * time.delta_t()
                                      : std::numeric_limits<amrex::Real>::max();
        const amrex::Real cfd_stop = amrex::min(stop1, stop2);
        const amrex::Real cfd_start = time.current_time();
        const amrex::Real cfd_sim = cfd_stop - cfd_start - 1.0e-6_rt;

        // Ensure that the user specified stop_time is not shorter than CFD sim
        const auto& tf = data.meta().ext_data;
        const amrex::Real fast_sim = tf.stop_time - tf.start_time;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fast_sim > cfd_sim,
            "OpenFAST simulation time is shorter than AMR-Wind duration");
    }
};

template <typename SrcTrait>
struct UpdatePosOp<TurbineFast, SrcTrait>
{
    void operator()(typename TurbineFast::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::UpdatePosOp<TurbineFast>");
        update_pos_op<TurbineFast::DataType>(data);
    }
};

template <typename SrcTrait>
struct UpdateVelOp<TurbineFast, SrcTrait>
{
    void operator()(typename TurbineFast::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::UpdateVelOp<TurbineFast>");
        update_vel_op<TurbineFast::DataType>(data);
    }
};

template <typename SrcTrait>
struct ComputeForceOp<TurbineFast, SrcTrait>
{
    void operator()(typename TurbineFast::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::ComputeForceOp<TurbineFast>");
        compute_force_op<TurbineFast::DataType>(data);
    }
};

template <typename SrcTrait>
struct ProcessOutputsOp<TurbineFast, SrcTrait>
{
private:
    typename TurbineFast::DataType& m_data;

    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir;

    //! NetCDF output filename for this turbine
    std::string m_nc_filename;

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename TurbineFast::DataType& data)
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

#endif /* TURBINE_FAST_OPS_H */
