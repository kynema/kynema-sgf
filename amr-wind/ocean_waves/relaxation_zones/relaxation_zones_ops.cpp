#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves::relaxation_zones {

void read_inputs(
    RelaxZonesBaseData& wdata,
    OceanWavesInfo& /*unused*/,
    const ::amr_wind::utils::MultiParser& pp)
{
    // Free surface zero level
    pp.query("zero_sea_level", wdata.zsl);
    pp.query("water_depth", wdata.water_depth);

    // Wave generation/absorption parameters
    pp.query("relax_zone_gen_length", wdata.gen_length);
    if (pp.contains("relax_zone_out_length")) {
        wdata.has_beach = false;
        wdata.init_wave_field = true;
        pp.query("relax_zone_out_length", wdata.beach_length);
    } else {
        pp.query("numerical_beach_length", wdata.beach_length);
        wdata.has_beach = true;
        pp.query("numerical_beach_length_factor", wdata.beach_length_factor);
        pp.query("initialize_wave_field", wdata.init_wave_field);
    }
    pp.query("relax_zone_length_y", wdata.zone_length_y);

    pp.query("current", wdata.current);

    wdata.has_ramp = pp.contains("timeramp_period");
    if (wdata.has_ramp) {
        pp.get("timeramp_period", wdata.ramp_period);
    }

    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.add("water_level", wdata.zsl);
}

void init_data_structures(RelaxZonesBaseData& /*unused*/) {}

void update_target_vof(CFDSim& sim)
{
    const int nlevels = sim.repo().num_active_levels();
    const auto& ow_levelset = sim.repo().get_field("ow_levelset");
    auto& ow_vof = sim.repo().get_field("ow_vof");
    const auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto target_phi = ow_levelset(lev).const_arrays();
        auto target_volfrac = ow_vof(lev).arrays();
        const amrex::Real eps = 2.0_rt * std::cbrt(dx[0] * dx[1] * dx[2]);
        amrex::ParallelFor(
            ow_vof(lev), amrex::IntVect(2),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                target_volfrac[nbx](i, j, k) =
                    multiphase::levelset_to_vof(i, j, k, eps, target_phi[nbx]);
            });
    }
    amrex::Gpu::streamSynchronize();
}

void modify_target_fields_for_beach(
    CFDSim& sim, const RelaxZonesBaseData& wdata)
{
    AMREX_ALWAYS_ASSERT(wdata.has_beach == true);
    const int nlevels = sim.repo().num_active_levels();
    auto& ow_levelset = sim.repo().get_field("ow_levelset");
    auto& ow_vel = sim.repo().get_field("ow_velocity");
    const auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();
        auto target_ls_arrs = ow_levelset(lev).arrays();
        auto target_vel_arrs = ow_vel(lev).arrays();

        const amrex::Real gen_length = wdata.gen_length;
        const amrex::Real beach_length = wdata.beach_length;
        const amrex::Real zsl = wdata.zsl;
        const amrex::Real current = wdata.current;

        amrex::ParallelFor(
            ow_vel(lev), amrex::IntVect(3),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // x is for combining wave profiles and does not need to exceed
                // domain limits
                const amrex::Real x = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(
                        problo[0] + ((i + 0.5_rt) * dx[0]), problo[0]),
                    probhi[0]);
                // z is for distance function and needs to exceed domain limits
                // so norms can be calculated when converted to vof
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                auto target_ls = target_ls_arrs[nbx];
                auto target_vel = target_vel_arrs[nbx];

                // Create wave vector for generation, numerical beach
                const utils::WaveVec wave_sol{
                    target_vel(i, j, k, 0), target_vel(i, j, k, 1),
                    target_vel(i, j, k, 2), target_ls(i, j, k)};
                const utils::WaveVec outlet{current, 0.0_rt, 0.0_rt, zsl - z};

                // Harmonize between inlet/bulk profile and outlet profile
                const auto target_profile = utils::harmonize_profiles_1d(
                    x, problo[0], gen_length, probhi[0], beach_length, wave_sol,
                    wave_sol, outlet);

                // Update target fields based on harmonization
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    target_vel(i, j, k, dir) = target_profile[dir];
                }
                target_ls(i, j, k) = target_profile[3];
            });
    }
    amrex::Gpu::streamSynchronize();
}

void apply_relaxation_zones(CFDSim& sim, const RelaxZonesBaseData& wdata)
{
    const int nlevels = sim.repo().num_active_levels();
    const auto& ow_vof = sim.repo().get_field("ow_vof");
    const auto& ow_vel = sim.repo().get_field("ow_velocity");
    const auto& geom = sim.mesh().Geom();

    const auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();
    constexpr amrex::Real vof_tiny = 1.0e-12_rt;

    // Get time
    const auto& time = sim.time().new_time();
    const amrex::Real rampf =
        (wdata.has_ramp) ? utils::ramp(time, wdata.ramp_period) : 1.0_rt;

    auto& vof = sim.repo().get_field("vof");
    auto& velocity = sim.repo().get_field("velocity");
    auto& density = sim.repo().get_field("density");

    amr_wind::IntField* terrain_blank_ptr{nullptr};
    amr_wind::IntField* terrain_drag_ptr{nullptr};
    const bool terrain_exists = sim.repo().int_field_exists("terrain_blank");
    // Get fields to prevent forcing in or near underwater terrain
    if (terrain_exists) {
        terrain_blank_ptr = &sim.repo().get_int_field("terrain_blank");
        terrain_drag_ptr = &sim.repo().get_int_field("terrain_drag");
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();
        auto vel_arrs = velocity(lev).arrays();
        auto rho_arrs = density(lev).arrays();
        auto volfrac_arrs = vof(lev).arrays();
        const auto target_volfrac_arrs = ow_vof(lev).const_arrays();
        const auto target_vel_arrs = ow_vel(lev).const_arrays();

        const auto terrain_blank_flags =
            terrain_exists ? (*terrain_blank_ptr)(lev).const_arrays()
                           : amrex::MultiArray4<int const>();
        const auto terrain_drag_flags =
            terrain_exists ? (*terrain_drag_ptr)(lev).const_arrays()
                           : amrex::MultiArray4<int const>();

        const amrex::Real gen_length = wdata.gen_length;
        const amrex::Real beach_length = wdata.beach_length;
        const amrex::Real beach_length_factor = wdata.beach_length_factor;
        const amrex::Real zone_length_y = wdata.zone_length_y;
        const bool has_zone_y = zone_length_y > constants::EPS;
        const bool has_beach = wdata.has_beach;

        amrex::ParallelFor(
            velocity(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(
                        problo[0] + ((i + 0.5_rt) * dx[0]), problo[0]),
                    probhi[0]);
                const amrex::Real y = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(
                        problo[1] + ((j + 0.5_rt) * dx[1]), problo[1]),
                    probhi[1]);

                auto vel = vel_arrs[nbx];
                auto rho = rho_arrs[nbx];
                auto volfrac = volfrac_arrs[nbx];
                const auto target_volfrac = target_volfrac_arrs[nbx];
                const auto target_vel = target_vel_arrs[nbx];

                // Skip if in or near terrain
                bool in_or_near_terrain{false};
                if (terrain_exists) {
                    in_or_near_terrain =
                        (terrain_blank_flags[nbx](i, j, k) == 1 ||
                         terrain_drag_flags[nbx](i, j, k) == 1);
                }

                // Get gamma for each possible direction
                const amrex::Real Gamma_xlo =
                    utils::gamma_generate(x - problo[0], gen_length);
                const amrex::Real Gamma_xhi = utils::gamma_absorb(
                    x - (probhi[0] - beach_length), beach_length,
                    beach_length_factor);
                amrex::Real Gamma_ylo = 1.0_rt;
                amrex::Real Gamma_yhi = 1.0_rt;
                if (has_zone_y) {
                    Gamma_ylo =
                        utils::gamma_generate(y - problo[1], zone_length_y);
                    Gamma_yhi = utils::gamma_absorb(
                        y - (probhi[1] - zone_length_y), zone_length_y, 1.0_rt);
                    const amrex::Real Gamma_y_to_xhi = utils::gamma_generate(
                        x - (probhi[0] - 1.5_rt * beach_length),
                        0.5_rt * beach_length);
                    if (has_beach) {
                        Gamma_ylo =
                            amrex::max<amrex::Real>(Gamma_ylo, Gamma_y_to_xhi);
                        Gamma_yhi =
                            amrex::max<amrex::Real>(Gamma_yhi, Gamma_y_to_xhi);
                    }
                }
                const amrex::Real Gamma = amrex::min<amrex::Real>(
                    amrex::min<amrex::Real>(Gamma_xhi, Gamma_xlo),
                    amrex::min<amrex::Real>(Gamma_yhi, Gamma_ylo));

                // Skip if Gamma is close enough to 1
                bool outside_zones = Gamma + constants::EPS >= 1.0_rt;

                if (!(outside_zones || in_or_near_terrain)) {
                    // Check for in beach, needed for velocity forcing
                    const bool in_beach =
                        has_beach && x + beach_length >= probhi[0];

                    // Target already harmonized: inlet/bulk with outlet profile
                    const utils::WaveVec target_profile{
                        target_vel(i, j, k, 0), target_vel(i, j, k, 1),
                        target_vel(i, j, k, 2), target_volfrac(i, j, k)};

                    // Nudge solution toward target
                    amrex::Real new_vof = utils::combine_linear(
                        Gamma, target_profile[3], volfrac(i, j, k));
                    new_vof = (new_vof > 1.0_rt - vof_tiny)
                                  ? 1.0_rt
                                  : (new_vof < vof_tiny ? 0.0_rt : new_vof);
                    const amrex::Real dvf = new_vof - volfrac(i, j, k);
                    volfrac(i, j, k) += rampf * dvf;
                    // Liquid velocity forced only where velocity is known
                    //  - in most of domain, that is where target vof is nonzero
                    //  - in numerical beach, that is anywhere
                    amrex::Real fvel_liq =
                        (target_volfrac(i, j, k) > vof_tiny) ? 1.0_rt : 0.0_rt;
                    fvel_liq = in_beach ? 1.0_rt : fvel_liq;
                    amrex::Real rho_ = (rho1 * volfrac(i, j, k)) +
                                       (rho2 * (1.0_rt - volfrac(i, j, k)));
                    for (int n = 0; n < vel.nComp(); ++n) {
                        // Get updated liquid velocity
                        amrex::Real vel_liq = vel(i, j, k, n);
                        const amrex::Real dvel_liq =
                            utils::combine_linear(
                                Gamma, target_profile[n], vel_liq) -
                            vel_liq;
                        vel_liq += rampf * fvel_liq * dvel_liq;
                        // If liquid was added, that liquid has target_vel
                        amrex::Real integrated_vel_liq =
                            volfrac(i, j, k) * vel_liq;
                        integrated_vel_liq +=
                            rampf * fvel_liq *
                            amrex::max<amrex::Real>(0.0_rt, dvf) *
                            (target_profile[n] - vel(i, j, k, n));
                        // Update overall velocity using momentum
                        vel(i, j, k, n) = (rho1 * integrated_vel_liq +
                                           rho2 * (1.0_rt - volfrac(i, j, k)) *
                                               vel(i, j, k, n)) /
                                          rho_;
                    }

                    // Make sure that density is updated before entering the
                    // solution
                    rho(i, j, k) = (rho1 * volfrac(i, j, k)) +
                                   (rho2 * (1.0_rt - volfrac(i, j, k)));
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    vof.fillpatch(time);
    velocity.fillpatch(time);
    density.fillpatch(time);
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const RelaxZonesBaseData& meta,
    const OceanWavesInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const RelaxZonesBaseData& meta,
    const OceanWavesInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace amr_wind::ocean_waves::relaxation_zones
