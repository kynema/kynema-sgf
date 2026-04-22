#include <algorithm>

#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"

namespace amr_wind::actuator::utils {

std::set<int> determine_influenced_procs(
    const amrex::AmrCore& mesh, const amrex::RealBox& rbx)
{
    std::set<int> procs;
    const int finest_level = mesh.finestLevel();
    const int nlevels = mesh.finestLevel() + 1;
    auto bx = amr_wind::utils::realbox_to_box(rbx, mesh.Geom(0));

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ba = mesh.boxArray(lev);
        const auto& dm = mesh.DistributionMap(lev);

        // Get all possible intersections at this level
        const auto& isects = ba.intersections(bx);

        // Extract the processor ranks
        for (const auto& is : isects) {
            procs.insert(dm[is.first]);
        }

        if (lev < finest_level) {
            bx = bx.refine(mesh.refRatio(lev));
        }
    }

    return procs;
}

void determine_root_proc(ActInfo& info, amrex::Vector<int>& act_proc_count)
{
    auto& plist = info.procs;
    bool assigned = false;

    // If any of the influenced procs is free (i.e., doesn't have a turbine
    // assigned to it) elect it as the root proc for this turbine and return
    // early.
    for (auto ip : plist) {
        if (act_proc_count[ip] < 1) {
            info.root_proc = ip;
            ++act_proc_count[ip];
            assigned = true;
            break;
        }
    }

    // If we found a root proc there is nothing more to do, so return early
    if (assigned) {
        const int iproc = amrex::ParallelDescriptor::MyProc();
        auto in_proc = info.procs.find(iproc);
        info.actuator_in_proc = (in_proc != info.procs.end());
        info.is_root_proc = (info.root_proc == iproc);

        // By default we request all processes where turbine is active to have
        // velocities sampled. Individual actuator instances can override this
        info.sample_vel_in_proc = info.actuator_in_proc;
        return;
    }

    // If we have reached here, then we have more turbines than processes
    // available. We will assign the current turbine to the process that is
    // managing the lowest number of turbines.

    // Determine the MPI rank that contains the fewest turbines
    auto it = std::ranges::min_element(act_proc_count);
    // Make it the root process for this turbine
    info.root_proc =
        static_cast<int>(std::distance(act_proc_count.begin(), it));
    // Make sure the root process is part of the process list
    plist.insert(info.root_proc);
    // Increment turbine count with the global tracking array
    ++(*it);

    {
        const int iproc = amrex::ParallelDescriptor::MyProc();
        auto in_proc = info.procs.find(iproc);
        info.actuator_in_proc = (in_proc != info.procs.end());
        info.is_root_proc = (info.root_proc == iproc);

        // By default we request all processes where turbine is active to have
        // velocities sampled. Individual actuator instances can override this
        info.sample_vel_in_proc = info.actuator_in_proc;
    }
}

} // namespace amr_wind::actuator::utils
