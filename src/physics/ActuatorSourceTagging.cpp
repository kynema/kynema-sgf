#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/ActuatorSourceTagging.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

ActuatorSourceTagging::ActuatorSourceTagging(CFDSim& sim) : m_repo(sim.repo())
{
    auto& pseqn = sim.pde_manager().register_transport_pde("PassiveScalar");
    m_tracer = &(pseqn.fields().field);

    amrex::ParmParse pp("ActuatorSourceTagging");
    pp.query("actuator_source_threshold", m_src_threshold);
}

void ActuatorSourceTagging::initialize_fields(
    int level, const amrex::Geometry& /*geom*/)
{
    (*m_tracer)(level).setVal(0.0_rt);
}

void ActuatorSourceTagging::post_init_actions()
{
    m_has_act_src = m_repo.field_exists("actuator_src_term");
    m_has_iblank = m_repo.field_exists("iblank_cell");

    if (m_has_act_src) {
        m_act_src = &(m_repo.get_field("actuator_src_term"));
        AMREX_ALWAYS_ASSERT(m_act_src->num_grow() <= m_tracer->num_grow());
    }

    if (m_has_iblank) {
        m_iblank = &(m_repo.get_int_field("iblank_cell"));
        AMREX_ALWAYS_ASSERT(m_iblank->num_grow() <= m_tracer->num_grow());
    }
}

void ActuatorSourceTagging::post_advance_work()
{
    if (!m_has_act_src && !m_has_iblank) {
        amrex::Print()
            << "Warning ActuatorSourceTagging activated but neither actuators "
               "or overset are being used"
            << '\n';
        return;
    }

    const amrex::Real src_threshold = m_src_threshold;
    for (int lev = 0; lev <= m_repo.mesh().finestLevel(); ++lev) {

        const auto& tracer_arrs = (*m_tracer)(lev).arrays();
        if (m_has_act_src) {
            const auto& src_arrs = (*m_act_src)(lev).const_arrays();
            amrex::ParallelFor(
                (*m_tracer)(lev), m_act_src->num_grow(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const auto src = src_arrs[nbx];
                    const amrex::Real srcmag = std::sqrt(
                        (src(i, j, k, 0) * src(i, j, k, 0)) +
                        (src(i, j, k, 1) * src(i, j, k, 1)) +
                        (src(i, j, k, 2) * src(i, j, k, 2)));

                    if (srcmag > src_threshold) {
                        tracer_arrs[nbx](i, j, k) = 1.0_rt;
                    }
                });
        }

        if (m_has_iblank) {
            const auto& iblank_arrs = (*m_iblank)(lev).const_arrays();
            const bool tag_fringe = m_tag_fringe;
            const bool tag_hole = m_tag_hole;
            amrex::ParallelFor(
                (*m_tracer)(lev), m_iblank->num_grow(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const auto ib = iblank_arrs[nbx](i, j, k);
                    if ((tag_fringe && (ib == -1)) || (tag_hole && (ib == 0))) {
                        tracer_arrs[nbx](i, j, k) = 1.0_rt;
                    }
                });
        }
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
