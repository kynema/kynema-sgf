#include "amr-wind/utilities/averaging/ReAveraging.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::averaging {

ReAveraging::ReAveraging(
    CFDSim& sim, const std::string& avgname, const std::string& fname)
    : m_field(sim.repo().get_field(fname))
    , m_average(sim.repo().declare_field(
          avg_name(m_field.name(), avgname),
          m_field.num_comp(),
          1, // 1 ghost cell to account for sampling
          1,
          m_field.field_location()))
{
    // Register default fillpatch operations
    m_average.set_default_fillpatch_bc(sim.time());
    // Do coarse/fine interpolations upon regrid
    m_average.fillpatch_on_regrid() = true;

    // Register average field with the IO manager
    auto& iomgr = sim.io_manager();
    iomgr.register_io_var(m_average.name());
}

const std::string& ReAveraging::average_field_name()
{
    return m_average.name();
}

void ReAveraging::operator()(
    const SimTime& time,
    const amrex::Real filter_width,
    const amrex::Real avg_time_interval,
    const amrex::Real elapsed_time)
{
    const amrex::Real filter =
        amrex::max(amrex::min(filter_width, elapsed_time), avg_time_interval);
    const amrex::Real factor =
        amrex::max<amrex::Real>(filter - avg_time_interval, 0.0_rt);

    const int ncomp = m_field.num_comp();
    const int nlevels = m_field.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ffab = m_field(lev);
        auto& afab = m_average(lev);

        const auto& fldarrs = ffab.const_arrays();
        const auto& avgarrs = afab.arrays();

        amrex::ParallelFor(
            ffab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                for (int n = 0; n < ncomp; ++n) {
                    const amrex::Real fval = fldarrs[nbx](i, j, k, n);
                    const amrex::Real aval = avgarrs[nbx](i, j, k, n);

                    avgarrs[nbx](i, j, k, n) =
                        (aval * factor + fval * avg_time_interval) / filter;
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    m_average.fillpatch(time.new_time());
}

} // namespace amr_wind::averaging
