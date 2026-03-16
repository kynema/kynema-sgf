#include "amr-wind/equation_systems/temperature/source_terms/PerturbationForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

PerturbationForcing::PerturbationForcing(const CFDSim& sim)
    : m_sim(sim), m_time(sim.time()), m_mesh(sim.mesh())
{
    amrex::ParmParse pp("PerturbationForcing");
    pp.queryarr("start", m_start, 0, AMREX_SPACEDIM);
    pp.queryarr("end", m_end, 0, AMREX_SPACEDIM);
    pp.query("time_steps", m_time_index);
    pp.query("start_level", m_start_level);
    pp.query("pert_amplitude", m_pert_amplitude);
}

PerturbationForcing::~PerturbationForcing() = default;

void PerturbationForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate */,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& index = m_time.time_index();
    if (index % m_time_index == 0 && index > 0) {
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> start{
            m_start[0], m_start[1], m_start[2]};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> end{
            m_end[0], m_end[1], m_end[2]};
        amrex::RealBox pert_box(start.data(), end.data());
        const bool has_terrain =
            this->m_sim.repo().int_field_exists("terrain_blank");
        const auto* m_terrain_blank =
            has_terrain ? &this->m_sim.repo().get_int_field("terrain_blank")
                        : nullptr;
        const auto& blank_arr = has_terrain
                                    ? (*m_terrain_blank)(lev).const_array(mfi)
                                    : amrex::Array4<int>();
        const auto* m_terrain_height =
            has_terrain ? &this->m_sim.repo().get_field("terrain_height")
                        : nullptr;
        const auto& height_arr = has_terrain
                                     ? (*m_terrain_height)(lev).const_array(mfi)
                                     : amrex::Array4<amrex::Real>();
        const auto& geom = m_mesh.Geom(lev);
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        const amrex::Real pert_amplitude = m_pert_amplitude;
        amrex::ParallelForRNG(
            bx, [=] AMREX_GPU_DEVICE(
                    int i, int j, int k, const amrex::RandomEngine& engine) {
                const amrex::Real height_arr_cell =
                    (has_terrain) ? height_arr(i, j, k, 0) : 0.0_rt;
                const amrex::Real blank_arr_cell =
                    (has_terrain) ? blank_arr(i, j, k, 0) : 0.0_rt;
                const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
                const amrex::Real z = amrex::max<amrex::Real>(
                    prob_lo[2] + ((k + 0.5_rt) * dx[2]) - height_arr_cell,
                    0.5_rt * dx[2]);
                const amrex::RealVect point{x, y, z};
                if (pert_box.contains(point)) {
                    const amrex::Real pert_cell =
                        amrex::RandomNormal(0, pert_amplitude, engine);
                    src_term(i, j, k, 0) +=
                        (1.0_rt - blank_arr_cell) * pert_cell;
                }
            });
    }
}

} // namespace amr_wind::pde::temperature
