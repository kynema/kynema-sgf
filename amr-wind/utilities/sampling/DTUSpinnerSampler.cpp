#include "amr-wind/utilities/sampling/DTUSpinnerSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_ParmParse.H"

#ifdef AMR_WIND_USE_OPENFAST
#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/turbine_fast_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
#endif

namespace amr_wind {
namespace sampling {

DTUSpinnerSampler::DTUSpinnerSampler(const CFDSim& sim) : LidarSampler(sim) {}

void DTUSpinnerSampler::initialize(const std::string& key)
{
    // Initialize the sampling time to be the same as simulation time
    m_time_sampling = m_sim.time().current_time();

    /*
     * Use this as a guide to implement the reading of the input
     * variables for the scanning pattern
     */

    // Read in new inputs specific to this class
    amrex::ParmParse pp(key);

    // Inner prism initial theta
    pp.get("inner_prism_theta0", m_InnerPrism.theta0);

    // Inner prism rotation rate
    pp.get("inner_prism_rotrate", m_InnerPrism.rot);

    // Inner prism azimuth angle
    pp.get("inner_prism_azimuth", m_InnerPrism.azimuth);

    // Outer prism initial theta
    pp.get("outer_prism_theta0", m_OuterPrism.theta0);

    // Outer prism rotation rate
    pp.get("outer_prism_rotrate", m_OuterPrism.rot);

    // Outer prism azimuth angle
    pp.get("outer_prism_azimuth", m_OuterPrism.azimuth);

    // This is the center of the lidar scan (x, y, z) [m]
    pp.getarr("lidar_center", m_lidar_center);
    AMREX_ALWAYS_ASSERT(
        static_cast<int>(m_lidar_center.size()) == AMREX_SPACEDIM);

    // Scan time
    pp.get("scan_time", m_scan_time);

    // Number of samples per scan
    pp.get("num_samples", m_num_samples);

    // Beam length
    pp.get("beam_length", m_beam_length);

    // Beam points
    pp.query("beam_points", m_beam_points);

    // Hub yaw logical flag
    pp.query("fixed_yaw", m_fixed_yaw);

    // Hub roll logical flag
    pp.query("fixed_roll", m_fixed_roll);

    // Hub tilt logical flag
    pp.query("fixed_tilt", m_fixed_tilt);

    // Spinner mode, hub or fixed
    pp.query("mode", m_spinner_mode);

    // For hub-mounted, this is turbine label
    pp.query("turbine", m_turbine_label);

    // Print hub fast turbine values to debug
    pp.query("hub_debug", m_hub_debug);

#ifdef AMR_WIND_USE_OPENFAST

    if (m_spinner_mode == "hub") {
        amrex::Print() << "Spinner Lidar will be attached to OpenFAST actuator"
                       << std::endl;
    }

#else

    AMREX_ALWAYS_ASSERT(m_spinner_mode == "fixed");

#endif // AMR_WIND_USE_OPENFAST

    update_sampling_locations();
    check_bounds();
}

vs::Vector DTUSpinnerSampler::adjust_lidar_pattern(
    vs::Vector beamPt, amrex::Real yaw, amrex::Real pitch, amrex::Real roll)
{

    const vs::Vector angles(roll, pitch, yaw);
    auto beamPt_transform = sampling_utils::rotation(angles, beamPt);
    return beamPt_transform;
}

vs::Vector DTUSpinnerSampler::generate_lidar_pattern(
    PrismParameters InnerPrism, PrismParameters OuterPrism, amrex::Real time)
{
    vs::Vector axis(1, 0, 0);
    vs::Vector ground(0, 0, 1);

    amrex::Real innerTheta = InnerPrism.theta0 + (InnerPrism.rot * time * 360);
    amrex::Real outerTheta = OuterPrism.theta0 + (OuterPrism.rot * time * 360);

    // NOLINTBEGIN(readability-suspicious-call-argument)
    auto reflection_1 = sampling_utils::rotate_euler_vec(
        axis, innerTheta,
        sampling_utils::rotate_euler_vec(
            ground, -((InnerPrism.azimuth / 2) + 90), axis));

    auto reflection_2 = sampling_utils::rotate_euler_vec(
        axis, outerTheta,
        sampling_utils::rotate_euler_vec(ground, OuterPrism.azimuth / 2, axis));
    // NOLINTEND(readability-suspicious-call-argument)

    return sampling_utils::reflect(
        reflection_2, sampling_utils::reflect(reflection_1, axis));
}

void DTUSpinnerSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    // Need a box that contains the fill val position so that it does
    // not get excluded from the write.
    const int lev = 0;
    const auto dom_hi = m_sim.mesh().Geom(lev).Domain().bigEnd();
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    const amrex::IntVect fill_val_iv(AMREX_D_DECL(
        static_cast<int>(amrex::Math::floor((m_fill_val - plo[0]) * dxinv[0])),
        static_cast<int>(amrex::Math::floor((m_fill_val - plo[1]) * dxinv[1])),
        static_cast<int>(
            amrex::Math::floor((m_fill_val - plo[2]) * dxinv[2]))));
    const amrex::Box box_containing_fill_val(fill_val_iv, dom_hi);
    sampling_locations(sample_locs, box_containing_fill_val);

    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void DTUSpinnerSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    const amrex::Real ndiv = amrex::max(m_beam_points - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;
    // Loop per subsampling
    for (int k = 0; k < m_ntotal; ++k) {

        int offset = k * AMREX_SPACEDIM;

        // Loop per spatial dimension
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            dx[d] = (m_end[d + offset] - m_start[d + offset]) / ndiv;
        }

        for (int i = 0; i < m_beam_points; ++i) {
            const amrex::RealVect loc = {AMREX_D_DECL(
                m_start[0 + offset] + (i * dx[0]),
                m_start[1 + offset] + (i * dx[1]),
                m_start[2 + offset] + (i * dx[2]))};
            if (utils::contains(box, loc, plo, dxinv)) {
                sample_locs.push_back(loc, i + (k * m_beam_points));
            }
        }
    }
}

#ifdef AMR_WIND_USE_OPENFAST
void DTUSpinnerSampler::bcast_turbine(
    amrex::Array<amrex::Real, 18>& turbine_pack, int root_proc)
{
    BL_PROFILE("amr-wind::Sampling::DTUSpinnerSampler::bcast_turbine");

    amrex::ParallelDescriptor::Bcast(
        turbine_pack.begin(), turbine_pack.size(), root_proc,
        amrex::ParallelDescriptor::Communicator());

    for (int i = 0; i < 9; i++) {
        m_current_hub_orient[i] = turbine_pack[i];
        if (i < 3) {
            m_current_hub_abs_pos[i] = static_cast<float>(turbine_pack[i + 9]);
            m_current_hub_rot_vel[i] = static_cast<float>(turbine_pack[i + 12]);
            m_turbine_base_pos[i] = static_cast<float>(turbine_pack[i + 15]);
        }
    }
}

void DTUSpinnerSampler::get_turbine_data(const std::string& turbine_label)
{

    BL_PROFILE("amr-wind::Sampling::DTUSpinnerSampler::get_turbine_data");

    // Use Physics Manager to get actuators
    const auto& all_actuators =
        m_sim.physics_manager().get<actuator::Actuator>();
    actuator::ActuatorModel& lidar_act =
        all_actuators.get_act_bylabel(turbine_label);

    if (m_sim.time().current_time() == m_sim.time().start_time()) {
        std::string actlabel = lidar_act.label();
        amrex::Print() << "Spinner Lidar Attached to actuator: " << actlabel
                       << std::endl;
    }

    // Test cast to Line or Disk
    auto* actline = dynamic_cast<
        actuator::ActModel<actuator::TurbineFast, actuator::ActSrcLine>*>(
        &lidar_act);

    auto* actdisk = dynamic_cast<
        actuator::ActModel<actuator::TurbineFast, actuator::ActSrcDisk>*>(
        &lidar_act);

    bool testline{(actline != nullptr)};
    bool testdisk{(actdisk != nullptr)};

    // Read and broadcast data
    if (testline) {
        const auto& actdata = actline->meta().ext_data;
        const auto& info = actline->info();

        // Create buffer object
        amrex::Array<amrex::Real, 18> turbine_pack = {};

        // Pack, broadcast, then unpack
        for (int i = 0; i < 9; i++) {
            turbine_pack[i] = actdata.hub_orient[i];
            if (i < 3) {
                turbine_pack[i + 9] = actdata.hub_abs_pos[i];
                turbine_pack[i + 12] = actdata.hub_rot_vel[i];
                turbine_pack[i + 15] = actdata.base_pos[i];
            }
        }

        bcast_turbine(turbine_pack, info.root_proc);
    } else if (testdisk) {
        const auto& actdata = actdisk->meta().ext_data;
        const auto& info = actdisk->info();

        // Create buffer object
        amrex::Array<amrex::Real, 18> turbine_pack = {};

        // Pack, broadcast, then unpack
        for (int i = 0; i < 9; i++) {
            turbine_pack[i] = actdata.hub_orient[i];
            if (i < 3) {
                turbine_pack[i + 9] = actdata.hub_abs_pos[i];
                turbine_pack[i + 12] = actdata.hub_rot_vel[i];
                turbine_pack[i + 15] = actdata.base_pos[i];
            }
        }

        bcast_turbine(turbine_pack, info.root_proc);
    } else {
        amrex::Abort("DTUSpinnerSampler: Problem finding actuator");
    }
}
#endif // AMR_WIND_USE_OPENFAST

bool DTUSpinnerSampler::update_sampling_locations()
{
    BL_PROFILE(
        "amr-wind::Sampling::DTUSpinnerSampler::update_sampling_locations");

#ifdef AMR_WIND_USE_OPENFAST
    if (m_spinner_mode == "hub") {
        get_turbine_data(m_turbine_label);

        m_hub_location = vs::Vector(
            m_current_hub_abs_pos[0] + m_turbine_base_pos[0],
            m_current_hub_abs_pos[1] + m_turbine_base_pos[1],
            m_current_hub_abs_pos[2] + m_turbine_base_pos[2]);

        // TODO: Do we need an offset from the hub location
        // to lidar start along shaft axis? Same for static angle misalignment?

        m_lidar_center[0] = m_hub_location[0];
        m_lidar_center[1] = m_hub_location[1];
        m_lidar_center[2] = m_hub_location[2];

        m_hub_tilt = std::atan2(
            -m_current_hub_orient[6],
            std::sqrt(
                std::pow(m_current_hub_orient[7], 2.0_rt) +
                std::pow(m_current_hub_orient[8], 2.0_rt)));
        m_hub_roll =
            std::atan2(m_current_hub_orient[7], m_current_hub_orient[8]);
        m_hub_yaw =
            std::atan2(m_current_hub_orient[3], m_current_hub_orient[0]);
    }
#endif

    // Sampling called in post_advance so time should be new_time
    // Not current_time()
    amrex::Real time = m_sim.time().new_time();
    amrex::Real start_time = m_sim.time().start_time();
    amrex::Real dt_sim = m_sim.time().delta_t();
    const amrex::Real dt_s = m_scan_time / m_num_samples;
    amrex::Real start_diff = std::abs(time - start_time);

    // Initialize the sampling time to the first time in the simulation
    if (start_diff < std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt &&
        m_update_count == 0) {
        m_time_sampling = time;
        m_hub_location_init = m_hub_location;
    }

    amrex::Real ts_diff = time - m_time_sampling;

    // Correction for time mismatch
    int time_corr = (ts_diff > dt_s) ? int(ts_diff / dt_s) : 0;

    m_ns = int(dt_sim / dt_s) + time_corr;
    m_ntotal = int(std::ceil(dt_sim / dt_s));

    int n_size = AMREX_SPACEDIM * (m_ns);
    int n_totalsize = AMREX_SPACEDIM * (m_ntotal);

    if (m_hub_debug) {
        amrex::Print() << "ts_diff: " << ts_diff << "\tm_ns: " << m_ns
                       << "\tSpin Time: " << m_time_sampling
                       << "\tAMR Time: " << time << "\tn_size: " << n_size
                       << '\n';
    }

    // Resize these variables so they can store all the locations
    if (m_start.size() < n_totalsize) {
        m_start.resize(n_totalsize);
    }
    if (m_end.size() < n_totalsize) {
        m_end.resize(n_totalsize);
    }

#ifdef AMR_WIND_USE_OPENFAST
    if (m_hub_debug) {
        amrex::Print() << "Turbine Hub Pos: " << m_current_hub_abs_pos[0] << " "
                       << m_current_hub_abs_pos[1] << " "
                       << m_current_hub_abs_pos[2] << " " << std::endl;
        amrex::Print() << "Turbine Rot Vel: " << m_current_hub_rot_vel[0] << " "
                       << m_current_hub_rot_vel[1] << " "
                       << m_current_hub_rot_vel[2] << " " << std::endl;
        amrex::Print() << "Turbine Orient: " << m_current_hub_orient[0] << " "
                       << m_current_hub_orient[1] << " "
                       << m_current_hub_orient[2] << " " << std::endl;
        amrex::Print() << "Turbine Orient: " << m_current_hub_orient[3] << " "
                       << m_current_hub_orient[4] << " "
                       << m_current_hub_orient[5] << " " << std::endl;
        amrex::Print() << "Turbine Orient: " << m_current_hub_orient[6] << " "
                       << m_current_hub_orient[7] << " "
                       << m_current_hub_orient[8] << " " << std::endl;
        amrex::Print() << "Yaw:Tilt:Roll: " << m_hub_yaw << " " << m_hub_tilt
                       << " " << m_hub_roll << std::endl;
        amrex::Print() << "Last Yaw:Tilt:Roll: " << m_last_hub_yaw << " "
                       << m_last_hub_tilt << " " << m_last_hub_roll
                       << std::endl;
    }
#endif

    if (m_update_count > 0) {
        // Loop per subsampling
        for (int k = 0; k < m_ntotal; k++) {

            int offset = k * AMREX_SPACEDIM;

            if (k < m_ns) {
                amrex::Real step =
                    (start_diff >
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt)
                        ? 1.0_rt * k
                        : 0.0_rt;
                amrex::Real srat =
                    (start_diff >
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt)
                        ? step / m_ns
                        : 0.0_rt;

                // Unit vector in the direction of the beam
                auto beam_vector = generate_lidar_pattern(
                    m_InnerPrism, m_OuterPrism, m_time_sampling);

                // Use fancy trick to fix interp at 180 to -180 transition of
                // atan2
                amrex::Real step_hub_yaw =
                    m_last_hub_yaw +
                    ((std::fmod(
                          std::fmod(m_hub_yaw - m_last_hub_yaw, m_twopi) +
                              m_threepi,
                          m_twopi) -
                      m_pi) *
                     srat);
                amrex::Real step_hub_tilt =
                    m_last_hub_tilt +
                    ((std::fmod(
                          std::fmod(m_hub_tilt - m_last_hub_tilt, m_twopi) +
                              m_threepi,
                          m_twopi) -
                      m_pi) *
                     srat);
                amrex::Real step_hub_roll =
                    m_last_hub_roll +
                    ((std::fmod(
                          std::fmod(m_hub_roll - m_last_hub_roll, m_twopi) +
                              m_threepi,
                          m_twopi) -
                      m_pi) *
                     srat);

                // Rotate beam unit vector
                beam_vector = adjust_lidar_pattern(
                    beam_vector, m_fixed_yaw + (step_hub_yaw * m_radtodeg),
                    m_fixed_tilt + (step_hub_tilt * m_radtodeg),
                    m_fixed_roll + (step_hub_roll * m_radtodeg));

                // Interpolate lidar center
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    m_step_lidar_center[d] =
                        ((step / m_ns) *
                         (m_lidar_center[d] - m_last_lidar_center[d])) +
                        m_last_lidar_center[d];
                }

                // Beam start and end points
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    m_start[d + offset] = m_step_lidar_center[d];
                    m_end[d + offset] = m_step_lidar_center[d] +
                                        (beam_vector[d] * m_beam_length);
                }

                m_time_sampling += dt_s;
            } else {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    m_start[d + offset] = m_fill_val;
                    m_end[d + offset] = m_fill_val;
                }
            }
        }
    }

    m_last_hub_yaw = m_hub_yaw;
    m_last_hub_tilt = m_hub_tilt;
    m_last_hub_roll = m_hub_roll;
    m_last_lidar_center = m_lidar_center;

    m_update_count = m_update_count + 1;

    return true;
}

#ifdef AMR_WIND_USE_NETCDF

void DTUSpinnerSampler::define_netcdf_metadata(
    const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
    grp.put_attr("time_table", m_time_table);
    grp.put_attr("azimuth_table", m_azimuth_table);
    grp.put_attr("elevation_table", m_elevation_table);

    std::vector<amrex::Real> scan_time{m_scan_time};
    std::vector<amrex::Real> num_samples{m_num_samples};
    std::vector<amrex::Real> beam_length{m_beam_length};
    std::vector<int> beam_points{m_beam_points};
    grp.put_attr("scan_time", scan_time);
    grp.put_attr("num_samples", num_samples);
    grp.put_attr("beam_length", beam_length);
    grp.put_attr("beam_points", beam_points);

    std::vector<amrex::Real> fixed_tilt{m_fixed_tilt};
    std::vector<amrex::Real> fixed_yaw{m_fixed_yaw};
    std::vector<amrex::Real> fixed_roll{m_fixed_roll};
    grp.put_attr("fixed_tilt", fixed_tilt);
    grp.put_attr("fixed_yaw", fixed_yaw);
    grp.put_attr("fixed_roll", fixed_roll);

    std::vector<amrex::Real> innerprism_theta0{m_InnerPrism.theta0};
    std::vector<amrex::Real> innerprism_rot{m_InnerPrism.rot};
    std::vector<amrex::Real> innerprism_azimuth{m_InnerPrism.azimuth};
    std::vector<amrex::Real> outerprism_theta0{m_OuterPrism.theta0};
    std::vector<amrex::Real> outerprism_rot{m_OuterPrism.rot};
    std::vector<amrex::Real> outerprism_azimuth{m_OuterPrism.azimuth};
    grp.put_attr("innerprism_theta0", innerprism_theta0);
    grp.put_attr("innerprism_rot", innerprism_rot);
    grp.put_attr("innerprism_azimuth", innerprism_azimuth);
    grp.put_attr("outerprism_theta0", outerprism_theta0);
    grp.put_attr("outerprism_rot", outerprism_rot);
    grp.put_attr("outerprism_azimuth", outerprism_azimuth);

    grp.put_attr("spinner_mode", m_spinner_mode);
    grp.put_attr("turbine", m_turbine_label);

    grp.def_dim("ndim", AMREX_SPACEDIM);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});

    grp.def_var("points_x", NC_DOUBLE, {"num_time_steps", "num_points"});
    grp.def_var("points_y", NC_DOUBLE, {"num_time_steps", "num_points"});
    grp.def_var("points_z", NC_DOUBLE, {"num_time_steps", "num_points"});

    grp.def_dim("nang", 3);
    grp.def_var("rotor_angles_rad", NC_DOUBLE, {"num_time_steps", "nang"});

    grp.def_var("rotor_hub_pos", NC_DOUBLE, {"num_time_steps", "ndim"});

    std::vector<amrex::Real> fillnan{std::nan("1")};
    grp.var("rotor_hub_pos").put_attr("_FillValue", fillnan);
    grp.var("rotor_angles_rad").put_attr("_FillValue", fillnan);
    grp.var("points").put_attr("_FillValue", fillnan);
}

void DTUSpinnerSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    std::vector<size_t> starti{nt, 0};
    std::vector<size_t> counti{1, 0};
    SampleLocType sample_locs;
    sampling_locations(sample_locs);

    auto xyz = grp.var("points");
    auto xp = grp.var("points_x");
    auto yp = grp.var("points_y");
    auto zp = grp.var("points_z");
    count[1] = num_points();
    counti[1] = num_points();

    const auto& locs = sample_locs.locations();
    xyz.put(locs[0].begin(), start, count);

    auto n_samples = m_beam_points * m_ntotal;

    std::vector<amrex::Real> xlocs(n_samples, 0.0_rt);
    std::vector<amrex::Real> ylocs(n_samples, 0.0_rt);
    std::vector<amrex::Real> zlocs(n_samples, 0.0_rt);

    // Loop per subsampling
    for (int k = 0; k < m_ntotal; ++k) {
        for (int i = 0; i < m_beam_points; ++i) {
            xlocs[i + k * m_beam_points] = locs[i + k * m_beam_points][0];
            ylocs[i + k * m_beam_points] = locs[i + k * m_beam_points][1];
            zlocs[i + k * m_beam_points] = locs[i + k * m_beam_points][2];
        }
    }

    xp.put(xlocs.data(), starti, counti);
    yp.put(ylocs.data(), starti, counti);
    zp.put(zlocs.data(), starti, counti);

    auto angs = grp.var("rotor_angles_rad");
    std::vector<size_t> acount{1, 3};
    amrex::Array<amrex::Real, 3> rangs = {m_hub_tilt, m_hub_roll, m_hub_yaw};
    angs.put(rangs.begin(), start, acount);

    auto hpos = grp.var("rotor_hub_pos");
    std::vector<size_t> pcount{1, AMREX_SPACEDIM};
    amrex::Array<amrex::Real, 3> rhpos = {
        m_lidar_center[0], m_lidar_center[1], m_lidar_center[2]};
    hpos.put(rhpos.begin(), start, pcount);
}

#else

void DTUSpinnerSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void DTUSpinnerSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}

#endif

} // namespace sampling

template struct ::amr_wind::sampling::SamplerBase::Register<
    ::amr_wind::sampling::DTUSpinnerSampler>;

} // namespace amr_wind
