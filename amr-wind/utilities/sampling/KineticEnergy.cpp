#include "amr-wind/utilities/sampling/KineticEnergy.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include <AMReX_MultiFabUtil.H>
#include <utility>
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::kinetic_energy {

KineticEnergy::KineticEnergy(CFDSim& sim, std::string label)
    : m_sim(sim)
    , m_label(std::move(label))
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{}

KineticEnergy::~KineticEnergy() = default;

void KineticEnergy::initialize()
{
    BL_PROFILE("amr-wind::KineticEnergy::initialize");
    amrex::ParmParse pp(m_label);
    populate_output_parameters(pp);
    prepare_ascii_file();
}

amrex::Real KineticEnergy::calculate_kinetic_energy()
{
    BL_PROFILE("amr-wind::KineticEnergy::calculate_kinetic_energy");

    // integrated total Kinetic Energy
    amrex::Real Kinetic_energy = 0.0_rt;

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

        Kinetic_energy += amrex::ReduceSum(
            m_density(lev), m_velocity(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& den_arr,
                amrex::Array4<amrex::Real const> const& vel_arr,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real Kinetic_Energy_Fab = 0.0_rt;

                amrex::Loop(bx, [=, &Kinetic_Energy_Fab](int i, int j, int k) {
                    Kinetic_Energy_Fab +=
                        cell_vol * mask_arr(i, j, k) * den_arr(i, j, k) *
                        (vel_arr(i, j, k, 0) * vel_arr(i, j, k, 0) +
                         vel_arr(i, j, k, 1) * vel_arr(i, j, k, 1) +
                         vel_arr(i, j, k, 2) * vel_arr(i, j, k, 2));
                });
                return Kinetic_Energy_Fab;
            });
    }

    // total volume of grid on level 0
    const amrex::Real total_vol = geom[0].ProbDomain().volume();

    Kinetic_energy *= 0.5_rt / total_vol;

    amrex::ParallelDescriptor::ReduceRealSum(Kinetic_energy);

    return Kinetic_energy;
}

void KineticEnergy::output_actions()
{
    BL_PROFILE("amr-wind::KineticEnergy::output_actions");

    m_total_kinetic_energy = calculate_kinetic_energy();

    write_ascii();
}

void KineticEnergy::prepare_ascii_file()
{
    BL_PROFILE("amr-wind::KineticEnergy::prepare_ascii_file");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());

    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        f << "time_step time kinetic_energy" << '\n';
        f.close();
    }
}

void KineticEnergy::write_ascii()
{
    BL_PROFILE("amr-wind::KineticEnergy::write_ascii");

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
        f << m_sim.time().time_index() << std::scientific
          << std::setprecision(m_precision) << std::setw(m_width)
          << m_sim.time().new_time();
        f << std::setw(m_width) << m_total_kinetic_energy;
        f << '\n';
        f.close();
    }
}

} // namespace amr_wind::kinetic_energy
