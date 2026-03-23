#include "amr-wind/utilities/averaging/ReynoldsStress.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::averaging {

ReynoldsStress::ReynoldsStress(
    CFDSim& sim, const std::string& avgname, const std::string& fname)
    : m_field(sim.repo().get_field("velocity"))
    , m_average(sim.repo().get_field("velocity_mean_" + avgname))
    , m_stress(sim.repo().declare_field(
          stress_name(m_field.name(), avgname),
          6, // number of components of the reynolds stress tensor
          1, // Ghost cells
          1,
          m_field.field_location()))
    , m_re_stress(sim.repo().declare_field(
          re_stress_name(m_field.name(), avgname),
          6, // number of components of the reynolds stress tensor
          1, // Ghost cells
          1,
          m_field.field_location()))
{
    if (fname != "velocity") {
        amrex::Abort("ReynoldsStress only implemented for velocity field");
    }

    // Register default fillpatch operations
    m_stress.set_default_fillpatch_bc(sim.time());
    m_re_stress.set_default_fillpatch_bc(sim.time());

    // Do coarse/fine interpolations upon regrid
    m_stress.fillpatch_on_regrid() = true;

    // Register average field with the IO manager
    auto& iomgr = sim.io_manager();
    iomgr.register_io_var(m_stress.name());
    iomgr.register_io_var(m_re_stress.name());
}

const std::string& ReynoldsStress::average_field_name()
{
    return m_re_stress.name();
}

void ReynoldsStress::operator()(
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
        const auto& afab = m_average(lev);
        auto& sfab = m_stress(lev);
        auto& rfab = m_re_stress(lev);

        const auto& fldarrs = ffab.const_arrays();
        const auto& avgarrs = afab.arrays();
        const auto& stressarrs = sfab.arrays();
        const auto& restressarrs = rfab.arrays();

        amrex::ParallelFor(
            ffab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // The tensor index
                int mn = 0;
                for (int n = 0; n < ncomp; ++n) {
                    for (int m = n; m < ncomp; ++m) {
                        // AB
                        const amrex::Real fval2 =
                            fldarrs[nbx](i, j, k, m) * fldarrs[nbx](i, j, k, n);
                        // <A><B>
                        const amrex::Real aval2 =
                            avgarrs[nbx](i, j, k, m) * avgarrs[nbx](i, j, k, n);
                        // The current value
                        const amrex::Real avg = stressarrs[nbx](i, j, k, mn);
                        // The stress <AB>
                        stressarrs[nbx](i, j, k, mn) =
                            (avg * factor + fval2 * avg_time_interval) / filter;
                        // The Reynolds stress <ab>
                        restressarrs[nbx](i, j, k, mn) =
                            stressarrs[nbx](i, j, k, mn) - aval2;
                        ++mn;
                    }
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    m_stress.fillpatch(time.new_time());
    m_re_stress.fillpatch(time.new_time());
}

} // namespace amr_wind::averaging
