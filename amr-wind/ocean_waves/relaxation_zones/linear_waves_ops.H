#ifndef LINEAR_WAVES_OPS_H
#define LINEAR_WAVES_OPS_H

#include "amr-wind/ocean_waves/relaxation_zones/LinearWaves.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/relaxation_zones/relaxation_zones_ops.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves::ops {

template <>
struct ReadInputsOp<LinearWaves>
{
    void operator()(
        LinearWaves::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        relaxation_zones::read_inputs(wdata, info, pp);
        // Get gravity, assume negative z
        amrex::Vector<amrex::Real> gravity{0.0_rt, 0.0_rt, -9.81_rt};
        amrex::ParmParse pp_incflo("incflo");
        pp_incflo.queryarr("gravity", gravity);
        wdata.g = -gravity[2];

        pp.get("wave_length", wdata.wave_length);
        pp.get("wave_height", wdata.wave_height);
        pp.query("wave_phase_offset_radians", wdata.wave_phase_offset);
        if (!pp.contains("wave_phase_offset_radians")) {
            pp.query("wave_phase_offset_degrees", wdata.wave_phase_offset);
            wdata.wave_phase_offset *=
                std::numbers::pi_v<amrex::Real> / 180.0_rt;
        } else if (pp.contains("wave_phase_offset_degrees")) {
            amrex::Abort(
                "ReadInputsOp<LinearWaves> : wave phase offset is specified in "
                "both radians and degrees. Please use only one.");
        }
    }
};

template <>
struct InitDataOp<LinearWaves>
{
    void operator()(
        LinearWaves::DataType& data,
        int level,
        const amrex::Geometry& geom,
        bool multiphase_mode)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();

        Field* levelset{nullptr};
        if (multiphase_mode) {
            levelset = &sim.repo().get_field("levelset");
        }
        // cppcheck-suppress constVariableReference
        auto& velocity = sim.repo().get_field("velocity");

        const auto& problo = geom.ProbLoArray();
        const auto& probhi = geom.ProbHiArray();
        const auto& dx = geom.CellSizeArray();

        const auto& vel = velocity(level).arrays();
        const auto& phi_arrs = multiphase_mode
                                   ? (*levelset)(level).arrays()
                                   : amrex::MultiArray4<amrex::Real>();

        const amrex::Real zero_sea_level = wdata.zsl;
        const amrex::Real gen_length = wdata.gen_length;
        const amrex::Real beach_length = wdata.beach_length;
        const amrex::Real g = wdata.g;
        const bool has_beach = wdata.has_beach && multiphase_mode;
        const bool init_wave_field = wdata.init_wave_field || !multiphase_mode;

        const amrex::Real wave_height = wdata.wave_height;
        const amrex::Real wave_length = wdata.wave_length;
        const amrex::Real phase_offset = wdata.wave_phase_offset;
        const amrex::Real water_depth = wdata.water_depth;
        const amrex::Real current = wdata.current;

        amrex::ParallelFor(
            velocity(level), amrex::IntVect(3),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                const amrex::Real wave_number =
                    2.0_rt * std::numbers::pi_v<amrex::Real> / wave_length;
                const amrex::Real omega = std::pow(
                    wave_number * g * std::tanh(wave_number * water_depth),
                    0.5_rt);
                const amrex::Real phase = (wave_number * x) - phase_offset;

                // Wave profile
                const amrex::Real eta_w =
                    (wave_height / 2.0_rt * std::cos(phase)) + zero_sea_level;
                const amrex::Real u_w =
                    omega * wave_height / 2.0_rt *
                    std::cosh(
                        wave_number * (z - zero_sea_level + water_depth)) /
                    std::sinh(wave_number * water_depth) * std::cos(phase);
                const amrex::Real v_w = 0.0_rt;
                const amrex::Real w_w =
                    omega * wave_height / 2.0_rt *
                    std::sinh(
                        wave_number * (z - zero_sea_level + water_depth)) /
                    std::sinh(wave_number * water_depth) * std::sin(phase);
                const utils::WaveVec wave_sol{u_w, v_w, w_w, eta_w};

                // Quiescent profile
                const utils::WaveVec quiescent{
                    0.0_rt, 0.0_rt, 0.0_rt, zero_sea_level};

                // Specify initial state for each region of domain
                const auto bulk = init_wave_field ? wave_sol : quiescent;
                const auto outlet = has_beach ? quiescent : wave_sol;

                const auto local_profile = utils::harmonize_profiles_1d(
                    x, problo[0], gen_length, probhi[0], beach_length, wave_sol,
                    bulk, outlet);

                const amrex::Real phi = local_profile[3] - z;
                const amrex::Real cell_length_2D =
                    std::sqrt((dx[0] * dx[0]) + (dx[2] * dx[2]));

                if (phi + cell_length_2D >= 0) {
                    vel[nbx](i, j, k, 0) = current + local_profile[0];
                    vel[nbx](i, j, k, 1) = local_profile[1];
                    vel[nbx](i, j, k, 2) = local_profile[2];
                }
                if (multiphase_mode) {
                    phi_arrs[nbx](i, j, k) = phi;
                }
            });
        amrex::Gpu::streamSynchronize();
    }
};

template <>
struct UpdateTargetFieldsOp<LinearWaves>
{
    void operator()(LinearWaves::DataType& data, const amrex::Real time)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();

        // cppcheck-suppress constVariableReference
        auto& ow_levelset = sim.repo().get_field("ow_levelset");
        // cppcheck-suppress constVariableReference
        auto& ow_velocity = sim.repo().get_field("ow_velocity");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            const auto& phi = ow_levelset(lev).arrays();
            const auto& vel = ow_velocity(lev).arrays();

            const amrex::Real wave_height = wdata.wave_height;
            const amrex::Real wave_length = wdata.wave_length;
            const amrex::Real phase_offset = wdata.wave_phase_offset;
            const amrex::Real water_depth = wdata.water_depth;
            const amrex::Real zero_sea_level = wdata.zsl;
            const amrex::Real g = wdata.g;
            const amrex::Real current = wdata.current;

            amrex::ParallelFor(
                ow_velocity(lev), amrex::IntVect(3),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const amrex::Real x = amrex::max<amrex::Real>(
                        problo[0], problo[0] + ((i + 0.5_rt) * dx[0]));
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    const amrex::Real wave_number =
                        2.0_rt * std::numbers::pi_v<amrex::Real> / wave_length;
                    const amrex::Real omega = std::pow(
                        wave_number * g * std::tanh(wave_number * water_depth),
                        0.5_rt);
                    const amrex::Real phase =
                        (wave_number * (x - current * time)) - (omega * time) -
                        phase_offset;

                    const amrex::Real eta =
                        (wave_height / 2.0_rt * std::cos(phase)) +
                        zero_sea_level;

                    phi[nbx](i, j, k) = eta - z;

                    const amrex::Real cell_length_2D =
                        std::sqrt((dx[0] * dx[0]) + (dx[2] * dx[2]));
                    if (phi[nbx](i, j, k) + cell_length_2D >= 0) {
                        vel[nbx](i, j, k, 0) =
                            (omega * wave_height / 2.0_rt *
                             std::cosh(
                                 wave_number *
                                 (z - zero_sea_level + water_depth)) /
                             std::sinh(wave_number * water_depth) *
                             std::cos(phase)) +
                            current;
                        vel[nbx](i, j, k, 1) = 0.0_rt;
                        vel[nbx](i, j, k, 2) =
                            omega * wave_height / 2.0_rt *
                            std::sinh(
                                wave_number *
                                (z - zero_sea_level + water_depth)) /
                            std::sinh(wave_number * water_depth) *
                            std::sin(phase);
                    } else {
                        vel[nbx](i, j, k, 0) = 0.0_rt;
                        vel[nbx](i, j, k, 1) = 0.0_rt;
                        vel[nbx](i, j, k, 2) = 0.0_rt;
                    }
                });
        }
        amrex::Gpu::streamSynchronize();
    }
};

} // namespace amr_wind::ocean_waves::ops

#endif /* LINEAR_WAVES_OPS_H */
