#ifndef ACTUATOR_OPSI_H
#define ACTUATOR_OPSI_H

#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/CFDSim.H"

#include <algorithm>

namespace amr_wind::actuator::ops {

template <typename T>
void determine_influenced_procs(typename T::DataType& data)
{
    auto& info = data.info();
    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    if (info.root_proc > -1) {
        // During regrid, the influenced processes might have changed and might
        // no longer include the root proc. We insert it back to ensure that it
        // is always present on the list.
        info.procs.insert(info.root_proc);
    }

    const int iproc = amrex::ParallelDescriptor::MyProc();
    auto in_proc = info.procs.find(iproc);
    info.actuator_in_proc = (in_proc != info.procs.end());
    info.sample_vel_in_proc = info.actuator_in_proc;
}

template <typename T>
void determine_root_proc(
    typename T::DataType& data, amrex::Vector<int>& act_proc_count)
{
    auto& info = data.info();

    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    utils::determine_root_proc(data.info(), act_proc_count);
}

} // namespace amr_wind::actuator::ops

#endif /* ACTUATOR_OPSI_H */
