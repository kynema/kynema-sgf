#include "amr-wind/utilities/sampling/WaveEnergy.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include <AMReX_MultiFabUtil.H>
#include <utility>
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::wave_energy {

WaveEnergy::WaveEnergy(CFDSim& sim, std::string label)
    : m_sim(sim)
    , m_label(std::move(label))
    , m_velocity(sim.repo().get_field("velocity"))
    , m_vof(sim.repo().get_field("vof"))
{}

WaveEnergy::~WaveEnergy() = default;

void WaveEnergy::initialize()
{
    BL_PROFILE("amr-wind::WaveEnergy::initialize");
    // Default water level is 0
    amrex::Real w_lev = 0.0_rt;
    {
        amrex::ParmParse pp(m_label);
        populate_output_parameters(pp);
        pp.query("water_level", w_lev);
    }
    // Get gravity and density
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
    // Calculate depth
    const auto& geom = m_sim.repo().mesh().Geom();
    amrex::Real depth = w_lev - geom[0].ProbLo()[2];
    // Calculate volume scaling and PE offset
    m_escl = (geom[0].ProbHi()[0] - geom[0].ProbLo()[0]) *
             (geom[0].ProbHi()[1] - geom[0].ProbLo()[1]) * depth;
    m_pe_off = -0.5_rt * m_gravity[2] * depth;

    prepare_ascii_file();
}

amrex::Real WaveEnergy::calculate_kinetic_energy()
{
    BL_PROFILE("amr-wind::WaveEnergy::calculate_kinetic_energy");

    // integrated total wave Energy
    amrex::Real wave_ke = 0.0_rt;

    const int finest_level = m_velocity.repo().num_active_levels() - 1;
    const auto& geom = m_velocity.repo().mesh().Geom();

    for (int lev = 0; lev <= finest_level; lev++) {

        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), m_sim.mesh().refRatio(lev), 1,
                0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        wave_ke += amrex::ReduceSum(
            m_vof(lev), m_velocity(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr,
                amrex::Array4<amrex::Real const> const& vel_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real Wave_Energy_Fab = 0.0_rt;

                amrex::Loop(bx, [=, &Wave_Energy_Fab](int i, int j, int k) {
                    Wave_Energy_Fab +=
                        cell_vol * mask_arr(i, j, k) * 0.5_rt *
                        vof_arr(i, j, k) *
                        (vel_arr(i, j, k, 0) * vel_arr(i, j, k, 0) +
                         vel_arr(i, j, k, 1) * vel_arr(i, j, k, 1) +
                         vel_arr(i, j, k, 2) * vel_arr(i, j, k, 2));
                });
                return Wave_Energy_Fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(wave_ke);

    return wave_ke;
}

amrex::Real WaveEnergy::calculate_potential_energy()
{
    BL_PROFILE("amr-wind::WaveEnergy::calculate_potential_energy");

    // integrated total wave Energy
    amrex::Real wave_pe = 0.0_rt;
    // gravity constant
    const amrex::Real g = -m_gravity[2];

    const int finest_level = m_velocity.repo().num_active_levels() - 1;
    const auto& geom = m_velocity.repo().mesh().Geom();

    for (int lev = 0; lev <= finest_level; lev++) {

        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                m_sim.mesh().boxArray(lev + 1), m_sim.mesh().refRatio(lev), 1,
                0);
        } else {
            level_mask.define(
                m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
                1, 0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        const amrex::Real dz = geom[lev].CellSize()[2];
        const amrex::Real probloz = geom[lev].ProbLo()[2];

        wave_pe += amrex::ReduceSum(
            m_vof(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real Wave_Energy_Fab = 0.0_rt;

                amrex::Loop(bx, [=, &Wave_Energy_Fab](int i, int j, int k) {
                    // Crude model of liquid height in multiphase cells
                    amrex::Real kk =
                        (vof_arr(i, j, k + 1) > vof_arr(i, j, k)) ? k + 1 : k;
                    amrex::Real dir =
                        (vof_arr(i, j, k + 1) > vof_arr(i, j, k)) ? -1 : 1;
                    const amrex::Real zl =
                        probloz + ((kk + dir * 0.5_rt * vof_arr(i, j, k)) * dz);
                    Wave_Energy_Fab += cell_vol * mask_arr(i, j, k) *
                                       vof_arr(i, j, k) * g * zl;
                });
                return Wave_Energy_Fab;
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(wave_pe);

    return wave_pe;
}

void WaveEnergy::output_actions()
{
    BL_PROFILE("amr-wind::WaveEnergy::output_actions");

    m_wave_kinetic_energy = calculate_kinetic_energy() / m_escl;
    m_wave_potential_energy =
        (calculate_potential_energy() / m_escl) + m_pe_off;

    write_ascii();
}

void WaveEnergy::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::WaveEnergy::prepare_ascii_file");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        f << "time_step time kinetic_energy potential_energy" << '\n';
        f.close();
    }
}

void WaveEnergy::write_ascii()
{
    BL_PROFILE("amr-wind::WaveEnergy::write_ascii");

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::scientific
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        f << std::setw(m_width) << m_wave_kinetic_energy;
        f << std::setw(m_width) << m_wave_potential_energy;
        f << '\n';
        f.close();
    }
}

} // namespace amr_wind::wave_energy
