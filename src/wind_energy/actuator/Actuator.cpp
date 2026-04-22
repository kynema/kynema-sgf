#include <algorithm>
#include <memory>
#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

Actuator::Actuator(CFDSim& sim)
    : m_sim(sim)
    , m_act_source(sim.repo().declare_field("actuator_src_term", 3, 1))
{}

Actuator::~Actuator() = default;

void Actuator::pre_init_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::pre_init_actions");
    amrex::ParmParse pp(identifier());

    amrex::Vector<std::string> labels;
    pp.getarr("labels", labels);
    ioutils::assert_with_message(
        ioutils::all_distinct(labels),
        "Duplicates in " + identifier() + ".labels");

    const int nturbines = static_cast<int>(labels.size());

    if (nturbines > 50) {
        amrex::Print()
            << "WARNING: There are many turbines in this case. Please ensure "
               "you have consolidated common turbine options under a common "
               "prefix. For example, the following: "
            << '\n';
        amrex::Print() << "  Actuator.Turb1.type            = UniformCtDisk"
                       << '\n';
        amrex::Print() << "  Actuator.Turb1.epsilon         = 5.0 5.0 5.0"
                       << '\n';
        amrex::Print() << "  Actuator.Turb2.type            = UniformCtDisk"
                       << '\n';
        amrex::Print() << "  Actuator.Turb2.epsilon         = 5.0 5.0 5.0"
                       << '\n';
        amrex::Print() << "becomes: " << '\n';
        amrex::Print() << "  Actuator.UniformCtDisk.epsilon = 5.0 5.0 5.0"
                       << '\n';
        amrex::Print() << "  Actuator.Turb1.type            = UniformCtDisk"
                       << '\n';
        amrex::Print() << "  Actuator.Turb2.type            = UniformCtDisk"
                       << '\n';
    }

    int cnt_turbfast = 0;
    for (int i = 0; i < nturbines; ++i) {
        const std::string& tname = labels[i];
        const std::string& prefix = identifier() + "." + tname;
        amrex::ParmParse pp1(prefix);

        std::string type;
        pp.query("type", type);
        pp1.query("type", type);
        AMREX_ALWAYS_ASSERT(!type.empty());
        cnt_turbfast +=
            ((type == "TurbineFastLine") || (type == "TurbineFastDisk")) ? 1
                                                                         : 0;

        auto obj = ActuatorModel::create(type, m_sim, tname, i);

        const std::string default_prefix = identifier() + "." + type;
        utils::ActParser inp(default_prefix, prefix);

        obj->read_inputs(inp);
        m_actuators.emplace_back(std::move(obj));
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        cnt_turbfast <= amrex::ParallelDescriptor::NProcs(),
        "Number of OpenFAST turbines must match or be lower than the number of "
        "ranks");
}

void Actuator::post_init_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_init_actions");

    amrex::Vector<int> act_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    for (auto& act : m_actuators) {
        act->determine_root_proc(act_proc_count);
    }

    {
        // Sanity check that we have processed the turbines correctly
        int nact =
            std::accumulate(act_proc_count.begin(), act_proc_count.end(), 0);
        AMREX_ALWAYS_ASSERT(num_actuators() == nact);
    }

    for (auto& act : m_actuators) {
        act->init_actuator_source();
    }

    // Ensure velocity fills ghost cells prior to sampling
    auto& vel = m_sim.repo().get_field("velocity");
    for (int lev = 0; lev < m_sim.repo().num_active_levels(); ++lev) {
        // Just need the interior velocity ghost cells updated
        vel(lev).FillBoundary(m_sim.mesh().Geom()[lev].periodicity());
    }

    setup_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
    prepare_outputs();
}

void Actuator::post_regrid_actions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_regrid_actions");

    const bool compute_root_proc =
        std::ranges::any_of(m_actuators, [](const auto& act) {
            return act->info().root_proc == -1;
        });
    if (compute_root_proc) {
        amrex::Vector<int> act_proc_count(
            amrex::ParallelDescriptor::NProcs(), 0);
        for (auto& act : m_actuators) {
            auto& info = act->info();
            info.procs =
                utils::determine_influenced_procs(m_sim.mesh(), info.bound_box);
            utils::determine_root_proc(act->info(), act_proc_count);
        }
    }

    for (auto& act : m_actuators) {
        act->determine_influenced_procs();
    }

    setup_container();
}

void Actuator::pre_advance_work()
{
    BL_PROFILE("amr-wind::actuator::Actuator::pre_advance_work");

    m_container->reset_container();
    update_positions();
    update_velocities();
    compute_forces();
    compute_source_term();
    communicate_turbine_io();
}

void Actuator::communicate_turbine_io()
{
#ifdef AMR_WIND_USE_HELICS
    BL_PROFILE("amr-wind::actuator::Actuator::communicate_turbine_io");

    if (!m_sim.helics().is_activated()) {
        return;
    }
    // send power and yaw from root actuator proc to io proc
    const int ptag = 0;
    const int ytag = 1;
    const size_t size = 1;
    for (auto& ac : m_actuators) {
        if (ac->info().is_root_proc) {
            amrex::ParallelDescriptor::Send(
                &m_sim.helics().m_turbine_power_to_controller[ac->info().id],
                size, amrex::ParallelDescriptor::IOProcessorNumber(), ptag);
            amrex::ParallelDescriptor::Send(
                &m_sim.helics()
                     .m_turbine_wind_direction_to_controller[ac->info().id],
                size, amrex::ParallelDescriptor::IOProcessorNumber(), ytag);
        }
        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::ParallelDescriptor::Recv(
                &m_sim.helics().m_turbine_power_to_controller[ac->info().id],
                size, ac->info().root_proc, ptag);
            amrex::ParallelDescriptor::Recv(
                &m_sim.helics()
                     .m_turbine_wind_direction_to_controller[ac->info().id],
                size, ac->info().root_proc, ytag);
        }
    }
#endif
}

/** Set up the container for sampling velocities
 *
 *  Allocates memory and initializes the particles corresponding to actuator
 *  nodes for all turbines that influence the current MPI rank. This method is
 *  invoked once during initialization and during regrid step.
 */
void Actuator::setup_container()
{
    BL_PROFILE("amr-wind::actuator::Actuator::setup_container");

    const int ntotal = num_actuators();
    const int nlocal = static_cast<int>(std::count_if(
        m_actuators.begin(), m_actuators.end(),
        [](const std::unique_ptr<ActuatorModel>& obj) {
            return obj->info().sample_vel_in_proc;
        }));

    m_container = std::make_unique<ActuatorContainer>(m_sim.mesh(), nlocal);

    auto& pinfo = m_container->m_data;
    for (int i = 0, il = 0; i < ntotal; ++i) {
        if (m_actuators[i]->info().sample_vel_in_proc) {
            pinfo.global_id[il] = i;
            pinfo.num_pts[il] = m_actuators[i]->num_velocity_points();
            ++il;
        }
    }

    m_container->initialize_container();
}

/** Update actuator positions and sample velocities at new locations.
 *
 *  This method loops over all the turbines local to this MPI rank and updates
 *  the position vectors. These new locations are provided to the sampling
 *  container that samples velocities at these new locations.
 *
 *  \sa Actuator::update_velocities
 */
void Actuator::update_positions()
{
    BL_PROFILE("amr-wind::actuator::Actuator::update_positions");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];
        auto vpos =
            ::amr_wind::utils::slice(pinfo.position, ic, pinfo.num_pts[i]);
        m_actuators[ig]->update_positions(vpos);
        ic += pinfo.num_pts[i];
    }
    m_container->update_positions();

    // Sample velocities at the new locations
    const auto& vel = m_sim.repo().get_field("velocity");
    const auto& density = m_sim.repo().get_field("density");
    m_container->sample_fields(vel, density);
}

/** Provide updated velocities from container to actuator instances
 *
 *  \sa Acuator::update_positions
 */
void Actuator::update_velocities()
{
    BL_PROFILE("amr-wind::actuator::Actuator::update_velocities");
    auto& pinfo = m_container->m_data;
    for (int i = 0, ic = 0; i < pinfo.num_objects; ++i) {
        const auto ig = pinfo.global_id[i];

        const auto vel =
            ::amr_wind::utils::slice(pinfo.velocity, ic, pinfo.num_pts[i]);

        const auto density =
            ::amr_wind::utils::slice(pinfo.density, ic, pinfo.num_pts[i]);

        m_actuators[ig]->update_fields(vel, density);
        ic += pinfo.num_pts[i];
    }
}

/** Helper method to compute forces on all actuator components
 */
void Actuator::compute_forces()
{
    BL_PROFILE("amr-wind::actuator::Actuator::compute_forces");
    for (auto& ac : m_actuators) {
        if (ac->info().actuator_in_proc) {
            ac->compute_forces();
        }
    }
}

void Actuator::compute_source_term()
{
    BL_PROFILE("amr-wind::actuator::Actuator::compute_source_term");
    m_act_source.setVal(0.0_rt);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& sfab = m_act_source(lev);
        const auto& geom = m_sim.mesh().Geom(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(sfab); mfi.isValid(); ++mfi) {
            for (auto& ac : m_actuators) {
                if (ac->info().actuator_in_proc) {
                    ac->compute_source_term(lev, mfi, geom);
                }
            }
        }

        // Ensure actuator src fills ghost cells prior to use (post-processing)
        // Just need the interior velocity ghost cells updated
        sfab.FillBoundary(geom.periodicity());
    }
}

void Actuator::prepare_outputs()
{
    BL_PROFILE("amr-wind::actuator::Actuator::prepare_outputs");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string out_dir_prefix = post_dir + "/actuator";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }
    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_actuators) {
        if (ac->info().root_proc == iproc) {
            ac->prepare_outputs(sname);
        }
    }
}

void Actuator::post_advance_work()
{
    BL_PROFILE("amr-wind::actuator::Actuator::post_advance_work");

    const int iproc = amrex::ParallelDescriptor::MyProc();
    for (auto& ac : m_actuators) {
        if (ac->info().root_proc == iproc) {
            ac->write_outputs();
        }
    }
}

ActuatorModel& Actuator::get_act_bylabel(const std::string& actlabel) const
{
    BL_PROFILE("amr-wind::actuator::Actuator::get_act_bylabel");
    int thisid = 0; // Default to first actuator
    for (const auto& act : m_actuators) {
        std::string thislabel = act->label();
        if (thislabel == actlabel) {
            thisid = act->id();
        }
    }

    return *m_actuators.at(thisid);
}

template <typename T>
T* Actuator::get_actuator(std::string& key) const
{
    BL_PROFILE("amr-wind::actuator::Actuator::get_actuator");
    for (const auto& act : m_actuators) {
        std::string thislabel = act->label();
        if (thislabel == key) {
            int thisid = act->id();
            T* converted = dynamic_cast<T*>(*m_actuators.at(thisid));
            return converted;
        }
    }

    amrex::Abort("Could not find actuator");
}

} // namespace amr_wind::actuator
