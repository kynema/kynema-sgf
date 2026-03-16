#ifndef W2A_WAVES_OPS_H
#define W2A_WAVES_OPS_H

#include "amr-wind/ocean_waves/relaxation_zones/W2AWaves.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

#ifdef AMR_WIND_USE_W2A
namespace {
int evaluate_read_resize(
    int& ntime,
    bool& read_flag,
    bool& resize_flag,
    amrex::Real& wtime,
    amrex::Real& t_last,
    const int new_ntime,
    const amrex::Real wtinit,
    const amrex::Real wdt,
    const amrex::Real time)
{
    // Flag to indicate that data must be read twice for interpolation
    int double_data = 0;
    // Check if time indicates reading must take place
    if (new_ntime != ntime) {
        // New data is needed - reading should happen
        read_flag = true;
        // If time index has proceeded more than 1 step
        if (new_ntime > ntime + 1) {
            // Double reading is necessary
            double_data = 1;
        }
        // Time index for reading
        ntime = new_ntime;
        // Sim time to go with recorded data
        wtime = new_ntime * wdt - wtinit;
        // If double reading is deemed necessary, check for convenience
        if (double_data == 1 && std::abs(wtime - time) <= 1.0e-10_rt) {
            // Reading can be done just once, w2a fields replace ow fields
            double_data = 2;
        }
    }
    // Check if reading must take place for other reasons
    if (t_last < -1.0e-10_rt) {
        // Signifies initialization from scratch without waves or a restart
        read_flag = true;
        // Resizing needs to happen for the first time
        resize_flag = true;

        // Confirm that new time is not coincident with modes time step
        if (std::abs(wtime - time) > 1.0e-10_rt) {
            // Data must be read before and after wtime
            double_data = 1;
        } else {
            // Data need only be read once, w2a fields will replace ow fields
            double_data = 2;
        }

    } else if (std::abs(t_last) < 1.0e-10_rt) {
        // Signifies initialization with waves
        read_flag = true;
        // Resizing needs to happen for the first time
        resize_flag = true;

        // levelset and velocity fields are up-to-date at t=0
        // interpolation is ready to go
    }
    // Record latest time as 'last' for next timestep
    t_last = time;
    // Return flag regarding double reading
    return double_data;
}

int update_offset_timestep(const int ntime, const int n0)
{
    // Offending timestep (goes too far): ntime + offset
    // Farthest back timesteps are permitted to go: n0
    // Subtract offset by offending timestep, add back lower limit
    // new offset = offset - (ntime + offset) + n0
    return (-ntime + n0);
}

void populate_fields_all_levels(
    amr_wind::ocean_waves::W2AWaves::MetaType& wdata,
    amrex::Vector<amrex::Geometry>& geom_all,
    amr_wind::Field& lvs_field,
    amr_wind::Field& vel_field,
    int ntime_off = 0)
{

    // Get data from modes
    bool no_EOF =
        wdata.is_ocean
            ? wdata.c_rmodes.get_data(
                  wdata.ntime + wdata.n_offset + ntime_off, wdata.c_mX,
                  wdata.c_mY, wdata.c_mZ, wdata.c_mFS)
            : wdata.r_rmodes.get_data(
                  wdata.ntime + wdata.n_offset + ntime_off, wdata.r_mX,
                  wdata.r_mY, wdata.r_mZ, wdata.r_mFS, wdata.r_mAdd);
    // Navigate when end of file is reached
    if (!no_EOF) {
        // End of file detected, reset reading
        wdata.n_offset = update_offset_timestep(wdata.ntime, wdata.n_winit);
        // Print warning to screen
        amrex::Print() << "WARNING (waves2amr_ops): end of mode data file "
                          "detected, resetting to beginning of mode data.\n";
        // Read data again, now from a valid timestep
        no_EOF = wdata.is_ocean
                     ? wdata.c_rmodes.get_data(
                           wdata.ntime + wdata.n_offset + ntime_off, wdata.c_mX,
                           wdata.c_mY, wdata.c_mZ, wdata.c_mFS)
                     : wdata.r_rmodes.get_data(
                           wdata.ntime + wdata.n_offset + ntime_off, wdata.r_mX,
                           wdata.r_mY, wdata.r_mZ, wdata.r_mFS, wdata.r_mAdd);
        // If no valid data is detected at this point, abort
        if (!no_EOF) {
            amrex::Abort(
                "waves2amr_ops: end of mode data file detected after "
                "resetting to beginning; please evaluate HOS_init_time "
                "or HOS_init_timestep and check the length of the mode "
                "file.");
        }
    }

    // Convert to spatial data in vectors
    if (wdata.is_ocean) {
        modes_hosgrid::copy_complex(
            wdata.n0, wdata.n1, wdata.c_mFS, wdata.c_eta_mptr);
        modes_hosgrid::populate_hos_eta(
            wdata.c_rmodes, wdata.plan_vector, wdata.c_eta_mptr,
            wdata.sp_eta_vec);
    } else {
        modes_hosgrid::copy_real(
            wdata.n0, wdata.n1, wdata.r_mFS, wdata.r_eta_mptr);
        modes_hosgrid::populate_hos_eta(
            wdata.r_rmodes, wdata.plan_vector, wdata.r_eta_mptr,
            wdata.sp_eta_vec);
    }

    for (int iht = 0; iht < wdata.indvec.size(); ++iht) {
        // Get sample height
        amrex::Real ht = wdata.hvec[wdata.indvec[iht]];
        // Sample velocity
        if (wdata.is_ocean) {
            modes_hosgrid::populate_hos_vel(
                wdata.c_rmodes, ht, wdata.zsl, wdata.c_mX, wdata.c_mY,
                wdata.c_mZ, wdata.plan_vector, wdata.c_u_mptr, wdata.c_v_mptr,
                wdata.c_w_mptr, wdata.sp_u_vec, wdata.sp_v_vec, wdata.sp_w_vec,
                iht * wdata.n0 * wdata.n1);
        } else {
            modes_hosgrid::populate_hos_vel(
                wdata.r_rmodes, ht, wdata.zsl, wdata.r_mX, wdata.r_mY,
                wdata.r_mZ, wdata.r_mAdd, wdata.plan_vector, wdata.r_u_mptr,
                wdata.r_v_mptr, wdata.r_w_mptr, wdata.au_mptr, wdata.av_mptr,
                wdata.aw_mptr, wdata.sp_u_vec, wdata.sp_v_vec, wdata.sp_w_vec,
                iht * wdata.n0 * wdata.n1);
        }
    }

    // Interpolate to fields (vector of MultiFabs)
    interp_to_mfab::interp_eta_to_levelset_field(
        wdata.n0_sp, wdata.n1_sp, wdata.dx0, wdata.dx1, wdata.xlo0, wdata.xlo1,
        wdata.zsl, wdata.is_ocean, wdata.sp_eta_vec, lvs_field.vec_ptrs(),
        geom_all);
    interp_to_mfab::interp_velocity_to_field(
        wdata.n0_sp, wdata.n1_sp, wdata.dx0, wdata.dx1, wdata.xlo0, wdata.xlo1,
        wdata.is_ocean, wdata.indvec, wdata.hvec, wdata.sp_u_vec,
        wdata.sp_v_vec, wdata.sp_w_vec, vel_field.vec_ptrs(), geom_all);
}

void time_interpolate_wave_fields(
    amr_wind::Field& lvs_field,
    amr_wind::Field& vel_field,
    const amr_wind::Field& future_lvs_field,
    const amr_wind::Field& future_vel_field,
    const amrex::Real t_old,
    const amrex::Real t_new,
    const amrex::Real t_future)
{
    const amrex::Real b =
        (t_new - t_old) /
        (t_future - t_old + std::numeric_limits<amrex::Real>::epsilon());
    const amrex::Real a = 1.0_rt - b;
    amr_wind::field_ops::lincomb(
        lvs_field, a, lvs_field, 0, b, future_lvs_field, 0, 0, 1,
        lvs_field.num_grow());
    amr_wind::field_ops::lincomb(
        vel_field, a, vel_field, 0, b, future_vel_field, 0, 0, AMREX_SPACEDIM,
        vel_field.num_grow());
}

} // namespace
#endif

namespace amr_wind::ocean_waves::ops {

template <>
struct ReadInputsOp<W2AWaves>
{
    void operator()(
        // cppcheck-suppress constParameterReference
        W2AWaves::DataType& data,
        const ::amr_wind::utils::MultiParser& pp)
    {
// Check for W2A initialization
#ifndef AMR_WIND_USE_W2A
        // Assert Waves2AMR must be used for initial condition file
        amrex::Abort(
            "ocean_waves/W2AWaves: AMR-Wind was not built with Waves2AMR "
            "support; associated wave data cannot be processed for relaxation "
            "zones.");

        amrex::ignore_unused(data, pp);
#else
        auto& wdata = data.meta();
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);

        if (wdata.current > constants::TIGHT_TOL) {
            amrex::Abort(
                "Current is specified as nonzero, but current is not yet "
                "implemented for W2A Waves.");
        }

        pp.get("HOS_modes_filename", wdata.modes_file);
        pp.query("HOS_simulation_is_ocean", wdata.is_ocean);
        pp.query("HOS_init_timestep", wdata.ntime);
        if (!pp.contains("HOS_init_timestep")) {
            pp.query("HOS_init_time", wdata.t_winit);
        }
        pp.query("HOS_domain_offset_x", wdata.xlo0);
        pp.query("HOS_domain_offset_y", wdata.xlo1);

        // Default fftw_plan is deterministic
        std::string fftw_planner_flag{"estimate"};
        pp.query("fftw_planner_flag", fftw_planner_flag);

        amrex::Vector<amrex::Real> prob_lo_input(AMREX_SPACEDIM);
        amrex::ParmParse pp_geom("geometry");
        pp_geom.getarr("prob_lo", prob_lo_input);

        // Read user inputs to height vector
        int nheights = 0;
        int nh_above = 1;
        amrex::Real dz0 = 0.0_rt;
        pp.get("number_interp_points_in_z", nheights);
        pp.get("interp_spacing_at_surface", dz0);
        pp.query("number_interp_above_surface", nh_above);

        // Initialize mode reader
        bool file_exists{false};
        if (wdata.is_ocean) {
            file_exists = wdata.c_rmodes.initialize(wdata.modes_file, true);
        } else {
            file_exists = wdata.r_rmodes.initialize(wdata.modes_file, false);
        }

        // Abort if file does not exist
        if (!file_exists) {
            amrex::Abort(
                "Waves2AMR ReadInputsOp: modes file requested does not exist");
        }

        // Get dt of HOS data
        wdata.dt_modes = wdata.is_ocean ? wdata.c_rmodes.get_dtout()
                                        : wdata.r_rmodes.get_dtout();

        // Get initial time and timestep synced
        if (wdata.t_winit > 0.0_rt) {
            // If initial time was specified
            // Get time index near requested time
            wdata.ntime =
                wdata.is_ocean
                    ? wdata.c_rmodes.time2step(wdata.t_winit, wdata.ntime)
                    : wdata.r_rmodes.time2step(wdata.t_winit, wdata.ntime);
            ;
            // Sync time to time index
            wdata.t_winit = wdata.dt_modes * wdata.ntime;
            // Save first timestep
            wdata.n_winit = wdata.ntime;
        } else {
            // If initial timestep is being used
            wdata.t_winit = wdata.dt_modes * wdata.ntime;
            // Save first timestep
            wdata.n_winit = wdata.ntime;
        }

        amrex::Print() << "OceanWaves: Initializing Waves2AMR profile at time "
                       << wdata.t_winit << ", step " << wdata.n_winit
                       << std::endl;

        // Initialize variables to store modes
        amrex::Real depth{0.0_rt};
        if (wdata.is_ocean) {
            int vsize = wdata.c_rmodes.get_vector_size();
            amrex::Real initval = 0.0_rt;
            wdata.c_mX.resize(vsize, initval);
            wdata.c_mY.resize(vsize, initval);
            wdata.c_mZ.resize(vsize, initval);
            wdata.c_mFS.resize(vsize, initval);

            // Get dimensions of data
            wdata.n0 = wdata.c_rmodes.get_first_fft_dimension();
            wdata.n1 = wdata.c_rmodes.get_second_fft_dimension();
            wdata.n0_sp = wdata.c_rmodes.get_first_spatial_dimension();
            wdata.n1_sp = wdata.c_rmodes.get_second_spatial_dimension();
            // Get resolution
            wdata.dx0 = wdata.c_rmodes.get_xlen() / wdata.n0_sp;
            wdata.dx1 = wdata.c_rmodes.get_ylen() / wdata.n1_sp;
            // Get depth
            depth = wdata.c_rmodes.get_depth();
            // Get dimensional length
            wdata.dimL = wdata.c_rmodes.get_L();
            // Get nominal last timestep of data
            wdata.n_wstop = static_cast<int>(
                (wdata.c_rmodes.get_Tstop() + 1.0e-8_rt) / wdata.dt_modes);

            amrex::Print()
                << "OceanWaves: Waves2AMR parameters from HOS-Ocean file: \n"
                << "    nx " << wdata.n0_sp << "  dx " << wdata.dx0 << "  Lx "
                << wdata.c_rmodes.get_xlen() << std::endl
                << "    ny " << wdata.n1_sp << "  dy " << wdata.dx1 << "  Ly "
                << wdata.c_rmodes.get_ylen() << std::endl
                << "    depth " << depth << std::endl
                << "    output time interval " << wdata.dt_modes
                << "  nominal last time " << wdata.n_wstop * wdata.dt_modes
                << std::endl;
        } else {
            int vsize = wdata.r_rmodes.get_vector_size();
            int vasize = wdata.r_rmodes.get_addl_vector_size();
            amrex::Real initval = 0.0_rt;
            wdata.r_mX.resize(vsize, initval);
            wdata.r_mY.resize(vsize, initval);
            wdata.r_mZ.resize(vsize, initval);
            wdata.r_mFS.resize(vsize, initval);
            wdata.r_mAdd.resize(vasize, initval);

            // Get dimensions of data
            wdata.n0 = wdata.r_rmodes.get_first_fft_dimension();
            wdata.n1 = wdata.r_rmodes.get_second_fft_dimension();
            wdata.n2 = wdata.r_rmodes.get_third_dimension();
            wdata.n0_sp = wdata.r_rmodes.get_first_spatial_dimension();
            wdata.n1_sp = wdata.r_rmodes.get_second_spatial_dimension();
            // Get resolution (slightly different for NWT)
            wdata.dx0 = wdata.r_rmodes.get_xlen() / (wdata.n0_sp - 1);
            wdata.dx1 =
                wdata.r_rmodes.get_ylen() / amrex::max(wdata.n1_sp - 1, 1);
            // Get depth
            depth = wdata.r_rmodes.get_depth();
            // Get dimensional length
            wdata.dimL = wdata.r_rmodes.get_L();
            // Get nominal last timestep of data
            wdata.n_wstop = static_cast<int>(
                (wdata.r_rmodes.get_Tstop() + 1.0e-8_rt) / wdata.dt_modes);

            amrex::Print()
                << "OceanWaves: Waves2AMR parameters from HOS-NWT file: \n"
                << "    nx " << wdata.n0_sp << "  dx " << wdata.dx0 << "  Lx "
                << wdata.r_rmodes.get_xlen() << std::endl
                << "    ny " << wdata.n1_sp << "  dy " << wdata.dx1 << "  Ly "
                << wdata.r_rmodes.get_ylen() << std::endl
                << "    nz " << wdata.n2 << "  depth " << depth << std::endl
                << "    output time interval " << wdata.dt_modes
                << "  nominal last time " << wdata.n_wstop * wdata.dt_modes
                << std::endl;
        }

        // Check if stop time is exceeded, introduce offset to ntime
        if (wdata.ntime + wdata.n_offset > wdata.n_wstop) {
            // If exceeding stop step, calculate new offset
            wdata.n_offset = update_offset_timestep(wdata.ntime, wdata.n_winit);
            // Print warning to screen
            amrex::Print()
                << "WARNING (waves2amr_ops): available mode data exceeded, "
                   "resetting to beginning of mode data.\n";
        }

        // Warning if depth does not correspond to simulation
        if (std::abs(depth - (wdata.zsl - prob_lo_input[2])) >
            1.0e-3_rt * depth) {
            amrex::Print()
                << "WARNING: Mismatch between water depths from AMR-Wind "
                   "domain and HOS data interpreted by Waves2AMR. \n   ^This "
                   "warning is not a concern when using waves as terrain.\n";
        }

        // Allocate pointers for FFTW
        if (wdata.is_ocean) {
            wdata.c_eta_mptr =
                modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
            wdata.c_u_mptr =
                modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
            wdata.c_v_mptr =
                modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
            wdata.c_w_mptr =
                modes_hosgrid::allocate_complex(wdata.n0, wdata.n1);
        } else {
            wdata.r_eta_mptr = modes_hosgrid::allocate_real(wdata.n0, wdata.n1);
            wdata.r_u_mptr = modes_hosgrid::allocate_real(wdata.n0, wdata.n1);
            wdata.r_v_mptr = modes_hosgrid::allocate_real(wdata.n0, wdata.n1);
            wdata.r_w_mptr = modes_hosgrid::allocate_real(wdata.n0, wdata.n1);
            wdata.au_mptr = modes_hosgrid::allocate_real(wdata.n0, 1);
            wdata.av_mptr = modes_hosgrid::allocate_real(wdata.n0, 1);
            wdata.aw_mptr = modes_hosgrid::allocate_real(wdata.n0, 1);
        }
        wdata.plan_out = modes_hosgrid::allocate_real(wdata.n0, wdata.n1);

        // Set up planner flag based on input
        auto plan_f = modes_hosgrid::planner_flags::estimate;
        if (fftw_planner_flag == "patient") {
            plan_f = modes_hosgrid::planner_flags::patient;
        } else if (fftw_planner_flag == "exhaustive") {
            plan_f = modes_hosgrid::planner_flags::exhaustive;
        } else if (fftw_planner_flag == "measure") {
            plan_f = modes_hosgrid::planner_flags::measure;
        } else if (!(fftw_planner_flag == "estimate")) {
            amrex::Print()
                << "WARNING (waves2amr_ops): invalid fftw_planner_flag "
                   "specified; defaulting to estimate (FFTW_ESTIMATE).\n";
        }
        // Set up plan for FFTW
        if (wdata.is_ocean) {
            wdata.plan_vector.emplace_back(
                modes_hosgrid::plan_ifftw(
                    wdata.n0, wdata.n1, wdata.c_eta_mptr, wdata.plan_out,
                    plan_f));
        } else {
            modes_hosgrid::plan_ifftw_nwt(
                wdata.n0, wdata.n1, wdata.plan_vector, wdata.r_eta_mptr,
                wdata.plan_out, plan_f);
        }

        // Create height vector for velocity mode conversion before
        // interpolation, with prob_lo as bottom
        int flag = interp_to_mfab::create_height_vector(
            wdata.hvec, nheights, dz0, wdata.zsl, prob_lo_input[2], nh_above);
        // Fail if flag indicates it should
        if (flag > 0) {
            amrex::Abort(
                "Waves2AMR ReadInputsOp: create_height_vector error, failure "
                "code " +
                std::to_string(flag));
        }

        amrex::Real z_double{0.0_rt};
        bool found_z_double{false};
        for (int n = 0; n < wdata.hvec.size() - 1; ++n) {
            const amrex::Real dz = wdata.hvec[n] - wdata.hvec[n + 1];
            if (!found_z_double && dz >= 2.0_rt * dz0) {
                z_double = 0.5_rt * (wdata.hvec[n] + wdata.hvec[n + 1]);
                found_z_double = true;
            }
        }

        amrex::Print()
            << "OceanWaves: Waves2AMR height vector for interpolation: \n"
            << "    start " << wdata.hvec[0] << "  end "
            << wdata.hvec[wdata.hvec.size() - 1] << "  num points "
            << wdata.hvec.size() << std::endl
            << "    initial dz " << wdata.hvec[0] - wdata.hvec[1]
            << "  first two grown dz "
            << wdata.hvec[nh_above] - wdata.hvec[nh_above + 1] << " "
            << wdata.hvec[nh_above + 1] - wdata.hvec[nh_above + 2]
            << "  final dz "
            << wdata.hvec[wdata.hvec.size() - 2] -
                   wdata.hvec[wdata.hvec.size() - 1]
            << std::endl;

        if (found_z_double) {
            amrex::Print()
                << "    interp spacing reaches double that of surface at "
                << z_double << std::endl;
        } else {
            amrex::Print() << "    interp spacing remains less than double "
                              "that of surface "
                           << std::endl;
        }

        // Always get modes on every processor for initialization
        bool no_EOF =
            wdata.is_ocean
                ? wdata.c_rmodes.get_data(
                      wdata.ntime + wdata.n_offset, wdata.c_mX, wdata.c_mY,
                      wdata.c_mZ, wdata.c_mFS)
                : wdata.r_rmodes.get_data(
                      wdata.ntime + wdata.n_offset, wdata.r_mX, wdata.r_mY,
                      wdata.r_mZ, wdata.r_mFS, wdata.r_mAdd);
        if (!no_EOF) {
            // End of file detected, reset reading
            wdata.n_offset = update_offset_timestep(wdata.ntime, wdata.n_winit);
            // Print warning to screen
            amrex::Print()
                << "WARNING (waves2amr_ops): end of mode data file "
                   "detected, resetting to beginning of mode data.\n";
            // Read data again, now from a valid timestep
            no_EOF =
                wdata.is_ocean
                    ? wdata.c_rmodes.get_data(
                          wdata.ntime + wdata.n_offset, wdata.c_mX, wdata.c_mY,
                          wdata.c_mZ, wdata.c_mFS)
                    : wdata.r_rmodes.get_data(
                          wdata.ntime + wdata.n_offset, wdata.r_mX, wdata.r_mY,
                          wdata.r_mZ, wdata.r_mFS, wdata.r_mAdd);
            // If no valid data is detected at this point, abort
            if (!no_EOF) {
                amrex::Abort(
                    "waves2amr_ops: end of mode data file detected after "
                    "resetting to beginning; please evaluate HOS_init_time "
                    "or HOS_init_timestep and check the length of the mode "
                    "file.");
            }
        }

        // Convert modes to spatial data
        if (wdata.is_ocean) {
            modes_hosgrid::copy_complex(
                wdata.n0, wdata.n1, wdata.c_mFS, wdata.c_eta_mptr);
            wdata.sp_eta_vec.resize(
                static_cast<size_t>(wdata.n0) * static_cast<size_t>(wdata.n1),
                0.0_rt);
            modes_hosgrid::populate_hos_eta(
                wdata.c_rmodes, wdata.plan_vector, wdata.c_eta_mptr,
                wdata.sp_eta_vec);
            // Mesh is not yet created, so get data at every height
            const auto n_hts = wdata.hvec.size();
            wdata.sp_u_vec.resize(
                static_cast<size_t>(wdata.n0 * wdata.n1) * n_hts);
            wdata.sp_v_vec.resize(
                static_cast<size_t>(wdata.n0 * wdata.n1) * n_hts);
            wdata.sp_w_vec.resize(
                static_cast<size_t>(wdata.n0 * wdata.n1) * n_hts);
            for (int iht = 0; iht < static_cast<int>(n_hts); ++iht) {
                // Get sample height
                amrex::Real ht = wdata.hvec[iht];
                // Sample velocity
                modes_hosgrid::populate_hos_vel(
                    wdata.c_rmodes, ht, wdata.zsl, wdata.c_mX, wdata.c_mY,
                    wdata.c_mZ, wdata.plan_vector, wdata.c_u_mptr,
                    wdata.c_v_mptr, wdata.c_w_mptr, wdata.sp_u_vec,
                    wdata.sp_v_vec, wdata.sp_w_vec, iht * wdata.n0 * wdata.n1);
            }
        } else {
            modes_hosgrid::copy_real(
                wdata.n0, wdata.n1, wdata.r_mFS, wdata.r_eta_mptr);
            wdata.sp_eta_vec.resize(
                static_cast<size_t>(wdata.n0) * static_cast<size_t>(wdata.n1),
                0.0_rt);
            modes_hosgrid::populate_hos_eta(
                wdata.r_rmodes, wdata.plan_vector, wdata.r_eta_mptr,
                wdata.sp_eta_vec);
            // Mesh is not yet created, so get data at every height
            const auto n_hts = wdata.hvec.size();
            wdata.sp_u_vec.resize(
                static_cast<size_t>(wdata.n0 * wdata.n1) * n_hts);
            wdata.sp_v_vec.resize(
                static_cast<size_t>(wdata.n0 * wdata.n1) * n_hts);
            wdata.sp_w_vec.resize(
                static_cast<size_t>(wdata.n0 * wdata.n1) * n_hts);
            for (int iht = 0; iht < static_cast<int>(n_hts); ++iht) {
                // Get sample height
                amrex::Real ht = wdata.hvec[iht];
                // Sample velocity
                modes_hosgrid::populate_hos_vel(
                    wdata.r_rmodes, ht, wdata.zsl, wdata.r_mX, wdata.r_mY,
                    wdata.r_mZ, wdata.r_mAdd, wdata.plan_vector, wdata.r_u_mptr,
                    wdata.r_v_mptr, wdata.r_w_mptr, wdata.au_mptr,
                    wdata.av_mptr, wdata.aw_mptr, wdata.sp_u_vec,
                    wdata.sp_v_vec, wdata.sp_w_vec, iht * wdata.n0 * wdata.n1);
            }
        }

        // Declare fields for HOS
        auto& w2a_levelset =
            data.sim().repo().declare_field("w2a_levelset", 1, 3, 2);
        auto& w2a_velocity = data.sim().repo().declare_field(
            "w2a_velocity", AMREX_SPACEDIM, 3, 2);

        // Extrapolation can work well when finer data is available
        w2a_levelset.set_default_fillpatch_bc(data.sim().time());
        w2a_velocity.set_default_fillpatch_bc(data.sim().time());
#endif
    }
}; // namespace ops

template <>
struct InitDataOp<W2AWaves>
{
    void
    // cppcheck-suppress constParameterReference
    operator()(
        W2AWaves::DataType& data,
        int level,
        const amrex::Geometry& geom,
        bool multiphase_mode)
    {

#ifdef AMR_WIND_USE_W2A
        auto& wdata = data.meta();

        auto& sim = data.sim();

        // Fill ow fields, then populate flow fields according to setup
        auto& ow_levelset = sim.repo().get_field("ow_levelset");
        auto& ow_velocity = sim.repo().get_field("ow_velocity");

        auto& velocity = sim.repo().get_field("velocity");
        // Set w2a fields to default values to prep for updates
        auto& w2a_levelset = sim.repo().get_field("w2a_levelset");
        auto& w2a_velocity = sim.repo().get_field("w2a_velocity");
        // Copy current fields to old w2a fields for time interp
        auto& w2a_lvs_old = w2a_levelset.state(amr_wind::FieldState::Old);
        auto& w2a_vel_old = w2a_velocity.state(amr_wind::FieldState::Old);

        Field* levelset{nullptr};
        if (multiphase_mode) {
            levelset = &sim.repo().get_field("levelset");
        }

        const auto& problo = geom.ProbLoArray();
        const auto& probhi = geom.ProbHiArray();
        const auto& dx = geom.CellSizeArray();

        // Set t_last to 0.0_rt to signify information read in
        wdata.t_last = 0.0_rt;

        // indvec is complete upon initialization (all heights every proc)
        amrex::Vector<int> indvec(wdata.hvec.size());
        std::iota(indvec.begin(), indvec.end(), 0);
        // Interpolate to MultiFabs (one level at a time)
        interp_to_mfab::interp_eta_to_levelset_multifab(
            wdata.n0_sp, wdata.n1_sp, wdata.dx0, wdata.dx1, wdata.xlo0,
            wdata.xlo1, wdata.zsl, wdata.is_ocean, wdata.sp_eta_vec,
            ow_levelset(level), problo, dx);
        interp_to_mfab::interp_velocity_to_multifab(
            wdata.n0_sp, wdata.n1_sp, wdata.dx0, wdata.dx1, wdata.xlo0,
            wdata.xlo1, wdata.is_ocean, indvec, wdata.hvec, wdata.sp_u_vec,
            wdata.sp_v_vec, wdata.sp_w_vec, ow_velocity(level), problo, dx);

        // Populate flow fields according to intended forcing and init setup
        const auto& ow_phi = ow_levelset(level).const_arrays();
        const auto& ow_vel = ow_velocity(level).const_arrays();
        const auto& vel = velocity(level).arrays();
        const auto& phi_arrs = multiphase_mode
                                   ? (*levelset)(level).arrays()
                                   : amrex::MultiArray4<amrex::Real>();

        const auto& w2a_phi = w2a_levelset(level).arrays();
        const auto& w2a_vel = w2a_velocity(level).arrays();
        const auto& w2a_phi_o = w2a_lvs_old(level).arrays();
        const auto& w2a_vel_o = w2a_vel_old(level).arrays();

        const amrex::Real gen_length = wdata.gen_length;
        const amrex::Real beach_length = wdata.beach_length;
        const amrex::Real y_length = wdata.zone_length_y;
        const amrex::Real zero_sea_level = wdata.zsl;

        const bool has_beach = wdata.has_beach && multiphase_mode;
        const bool init_wave_field = wdata.init_wave_field || !multiphase_mode;

        amrex::MultiFab phi_tmp_fab(
            velocity(level).boxArray(), velocity(level).DistributionMap(), 1,
            3);
        amrex::MultiFab vel_tmp_fab(
            velocity(level).boxArray(), velocity(level).DistributionMap(),
            AMREX_SPACEDIM, 3);

        const auto& phi_tmp = phi_tmp_fab.arrays();
        const auto& vel_tmp = vel_tmp_fab.arrays();

        amrex::ParallelFor(
            velocity(level), amrex::IntVect(3),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + 0.5_rt) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5_rt) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5_rt) * dx[2];

                // Wave profile
                const utils::WaveVec wave_sol{
                    ow_vel[nbx](i, j, k, 0), ow_vel[nbx](i, j, k, 1),
                    ow_vel[nbx](i, j, k, 2), ow_phi[nbx](i, j, k) + z};
                // Quiescent profile
                const utils::WaveVec quiescent{
                    0.0_rt, 0.0_rt, 0.0_rt, zero_sea_level};

                // Specify initial state for each region of domain
                const auto bulk = init_wave_field ? wave_sol : quiescent;
                const auto outlet = has_beach ? quiescent : wave_sol;

                const auto local_profile = utils::harmonize_profiles_2d(
                    x, y, problo[0], gen_length, probhi[0], beach_length,
                    problo[1], probhi[1], y_length, wave_sol, bulk, outlet);

                phi_tmp[nbx](i, j, k) = local_profile[3] - z;
                vel_tmp[nbx](i, j, k, 0) = local_profile[0];
                vel_tmp[nbx](i, j, k, 1) = local_profile[1];
                vel_tmp[nbx](i, j, k, 2) = local_profile[2];

                if (multiphase_mode) {
                    phi_arrs[nbx](i, j, k) = phi_tmp[nbx](i, j, k);
                }

                // Start old values at w2a solution (for the sake of time
                // interpolation)
                w2a_phi_o[nbx](i, j, k) = wave_sol[3] - z;
                w2a_vel_o[nbx](i, j, k, 0) = wave_sol[0];
                w2a_vel_o[nbx](i, j, k, 1) = wave_sol[1];
                w2a_vel_o[nbx](i, j, k, 2) = wave_sol[2];

                // Default w2a values matter for where no updates happen
                w2a_phi[nbx](i, j, k) = quiescent[3] - z;
                w2a_vel[nbx](i, j, k, 0) = quiescent[0];
                w2a_vel[nbx](i, j, k, 1) = quiescent[1];
                w2a_vel[nbx](i, j, k, 2) = quiescent[2];
            });
        amrex::Gpu::streamSynchronize();

        amrex::ParallelFor(
            velocity(level), amrex::IntVect(2),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real eps =
                    2.0_rt * std::cbrt(dx[0] * dx[1] * dx[2]);
                const amrex::Real vof_local =
                    multiphase::levelset_to_vof(i, j, k, eps, phi_tmp[nbx]);

                if (vof_local > 1.0_rt - constants::TIGHT_TOL) {
                    // Entire cell is liquid
                    vel[nbx](i, j, k, 0) = vel_tmp[nbx](i, j, k, 0);
                    vel[nbx](i, j, k, 1) = vel_tmp[nbx](i, j, k, 1);
                    vel[nbx](i, j, k, 2) = vel_tmp[nbx](i, j, k, 2);
                } else if (vof_local > 1.0e-3_rt) {
                    // Velocity of liquid becomes velocity of entire cell
                    // Interpolate using cell-centered values to liquid centroid
                    // Consider only z-direction to find centroid
                    const amrex::Real z_cell = problo[2] + (k + 0.5_rt) * dx[2];
                    const amrex::Real z_cell_below = z_cell - dx[2];
                    const amrex::Real z_cell_bottom = z_cell - 0.5_rt * dx[2];
                    const amrex::Real z_liq_center =
                        z_cell_bottom + 0.5_rt * vof_local * dx[2];
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        vel[nbx](i, j, k, n) =
                            vel_tmp[nbx](i, j, k - 1, n) +
                            (vel_tmp[nbx](i, j, k, n) -
                             vel_tmp[nbx](i, j, k - 1, n)) *
                                ((z_liq_center - z_cell_below) /
                                 (z_cell - z_cell_below + constants::EPS));
                    }
                }
            });
        amrex::Gpu::streamSynchronize();
#else
        amrex::ignore_unused(data, level, geom, multiphase_mode);
#endif
    }
}; // namespace ocean_waves

template <>
struct UpdateTargetFieldsOp<W2AWaves>
{
    // cppcheck-suppress constParameterReference
    void operator()(W2AWaves::DataType& data, const amrex::Real time)
    {

#ifdef AMR_WIND_USE_W2A
        auto& wdata = data.meta();
        auto& sim = data.sim();

        // Update ow fields every time
        auto& ow_levelset = sim.repo().get_field("ow_levelset");
        auto& ow_velocity = sim.repo().get_field("ow_velocity");
        // Update HOS fields when necessary
        auto& w2a_levelset = sim.repo().get_field("w2a_levelset");
        auto& w2a_velocity = sim.repo().get_field("w2a_velocity");
        // Old state of HOS fields used for time interpolation
        auto& w2a_lvs_old = w2a_levelset.state(amr_wind::FieldState::Old);
        auto& w2a_vel_old = w2a_velocity.state(amr_wind::FieldState::Old);

        // Proxy for multiphase mode of ocean waves
        bool vof_exists = sim.repo().field_exists("vof");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        // Get value for time interpolation
        amrex::Real t_last = wdata.t_last;

        // Check if new HOS data needs to be read
        bool read_flag = false;
        // Check if time indicates reading must take place
        int new_ntime =
            wdata.is_ocean
                ? wdata.c_rmodes.time2step(time + wdata.t_winit, wdata.ntime)
                : wdata.r_rmodes.time2step(time + wdata.t_winit, wdata.ntime);
        int double_data = evaluate_read_resize(
            wdata.ntime, read_flag, wdata.resize_flag, wdata.t, wdata.t_last,
            new_ntime, wdata.t_winit, wdata.dt_modes, time);
        // Check if stop time is exceeded, introduce offset to ntime
        if (read_flag) {
            // Need to only check when reading is occurring
            if (wdata.ntime + wdata.n_offset > wdata.n_wstop) {
                // If exceeding stop step, calculate new offset
                wdata.n_offset =
                    update_offset_timestep(wdata.ntime, wdata.n_winit);
                // Print warning to screen
                amrex::Print()
                    << "WARNING (waves2amr_ops): available mode data exceeded, "
                       "resetting to beginning of mode data.\n";
            }
        }
        // Resizing (assuming reading is taking place) must happen after regrid
        if (wdata.regrid_occurred) {
            // resize_flag remains true until resizing occurs, but
            // regrid_occurred resets every timestep
            wdata.resize_flag = true;
        }

        // Read HOS data if necessary based on time
        if (read_flag) {

            amrex::Print() << "OceanWaves: Waves2AMR reading from file: step "
                           << wdata.ntime << std::endl
                           << "    last interp time " << t_last << std::endl
                           << "    simulation time  " << time << std::endl
                           << "    new data time    " << wdata.t << std::endl
                           << "    double data read " << double_data
                           << std::endl;

            if (wdata.resize_flag) {
                // Reset flag
                wdata.resize_flag = false;
                // Flags for indicating overlap, assume none at first
                bool flag_z = false;
                bool flag_xlo = false;
                bool flag_xhi = false;
                bool flag_ylo = false;
                bool flag_yhi = false;
                // Get heights for this processor, check overlap in z
                flag_z =
                    (interp_to_mfab::get_local_height_indices(
                         wdata.indvec, wdata.hvec, ow_velocity.vec_ptrs(),
                         geom) == 1);
                // No overlap from heights definitely means no interp

                // Check lateral bounds (in x)
                flag_xlo =
                    (interp_to_mfab::check_lateral_overlap_lo(
                         wdata.gen_length, 0, ow_velocity.vec_ptrs(), geom) ==
                     1);
                // No overlap with gen region means no interp, unless ...
                if (!wdata.has_beach) {
                    // ... if overlap exists here, needing interp
                    flag_xhi =
                        (interp_to_mfab::check_lateral_overlap_hi(
                             wdata.beach_length, 0, ow_velocity.vec_ptrs(),
                             geom) == 1);
                }
                // Then in y, if y zones are present
                if (wdata.zone_length_y > constants::EPS) {
                    flag_ylo =
                        (interp_to_mfab::check_lateral_overlap_lo(
                             wdata.zone_length_y, 1, ow_velocity.vec_ptrs(),
                             geom) == 1);
                    flag_yhi =
                        (interp_to_mfab::check_lateral_overlap_hi(
                             wdata.zone_length_y, 1, ow_velocity.vec_ptrs(),
                             geom) == 1);
                }

                // Wave generation zones are irrelevant for single-phase mode,
                // indicated by nonexistence of vof
                if (!vof_exists) {
                    flag_xlo = true;
                    flag_xhi = true;
                }

                if (flag_z && (flag_xlo || flag_xhi || flag_ylo || flag_yhi)) {
                    // Interpolation is needed
                    wdata.do_interp = true;
                    // Do resizing
                    wdata.sp_eta_vec.resize(
                        static_cast<size_t>(wdata.n0) *
                            static_cast<size_t>(wdata.n1),
                        0.0_rt);
                    wdata.sp_u_vec.resize(
                        static_cast<size_t>(wdata.n0 * wdata.n1) *
                        wdata.indvec.size());
                    wdata.sp_v_vec.resize(
                        static_cast<size_t>(wdata.n0 * wdata.n1) *
                        wdata.indvec.size());
                    wdata.sp_w_vec.resize(
                        static_cast<size_t>(wdata.n0 * wdata.n1) *
                        wdata.indvec.size());
                    // Sizes will remain constant and need for interpolation
                    // will remain until a regrid occurs
                } else {
                    // No overlapping with spatial data or no overlapping with
                    // relaxation zones, interpolation can be skipped
                    wdata.do_interp = false;
                }
            }
            // Only perform reading where needed, communicate offset though
            amrex::ParallelDescriptor::ReduceIntMax(wdata.n_offset);

            // Previous W2A data is more recent than last interpolated data
            t_last = (wdata.ntime - 1) * wdata.dt_modes;

            if (double_data == 0) {
                // Previous W2A data needs to be copied into old fields
                amr_wind::field_ops::copy(
                    w2a_lvs_old, w2a_levelset, 0, 0, 1,
                    w2a_levelset.num_grow());
                amr_wind::field_ops::copy(
                    w2a_vel_old, w2a_velocity, 0, 0, AMREX_SPACEDIM,
                    w2a_velocity.num_grow());
            } else if (double_data == 1) {
                // Previous W2A data is unknown, read into old fields
                if (wdata.do_interp) {
                    populate_fields_all_levels(
                        wdata, geom, w2a_lvs_old, w2a_vel_old, -1);
                }

                // Average down to get fine information on coarse grid where
                // possible (may be unnecessary)
                for (int lev = nlevels - 1; lev > 0; --lev) {
                    amrex::average_down(
                        w2a_vel_old(lev), w2a_vel_old(lev - 1), 0,
                        AMREX_SPACEDIM, sim.mesh().refRatio(lev - 1));
                    amrex::average_down(
                        w2a_lvs_old(lev), w2a_lvs_old(lev - 1), 0, 1,
                        sim.mesh().refRatio(lev - 1));
                }
                // Fill patch to get correct ghost cells after average down
                w2a_vel_old.fillpatch(sim.time().new_time());
                w2a_lvs_old.fillpatch(sim.time().new_time());
            }
            // ow fields cancel when double_data == 2, no modification

            // After possible prior read, now read data for this ntime
            if (wdata.do_interp) {
                populate_fields_all_levels(
                    wdata, geom, w2a_levelset, w2a_velocity);
            }

            // Average down to get fine information on coarse grid where
            // possible and update ghost cells (may be unnecessary)
            for (int lev = nlevels - 1; lev > 0; --lev) {
                amrex::average_down(
                    w2a_velocity(lev), w2a_velocity(lev - 1), 0, AMREX_SPACEDIM,
                    sim.mesh().refRatio(lev - 1));
                w2a_velocity(lev - 1).FillBoundary(geom[lev - 1].periodicity());
                amrex::average_down(
                    w2a_levelset(lev), w2a_levelset(lev - 1), 0, 1,
                    sim.mesh().refRatio(lev - 1));
                w2a_levelset(lev - 1).FillBoundary(geom[lev - 1].periodicity());
            }
        }

        // Temporally interpolate at every timestep to get target solution
        time_interpolate_wave_fields(
            w2a_lvs_old, w2a_vel_old, w2a_levelset, w2a_velocity, t_last, time,
            wdata.t);
        // Copy to ow fields, which may be modified for beach
        field_ops::copy(
            ow_levelset, w2a_lvs_old, 0, 0, 1, ow_levelset.num_grow());
        field_ops::copy(
            ow_velocity, w2a_vel_old, 0, 0, AMREX_SPACEDIM,
            ow_velocity.num_grow());
#else
        amrex::ignore_unused(data, time);
#endif
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* WAVES2AMR_OPS_H */
