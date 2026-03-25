#ifndef Temperature_OPS_H
#define Temperature_OPS_H

#include "amr-wind/equation_systems/temperature/temperature.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind::pde {

/** Implementation of post-processing actions for temperature equation
 *  \ingroup tke_eqn
 */
template <>
struct PostSolveOp<Temperature>
{
    PostSolveOp(CFDSim& /*unused*/, PDEFields& fields) : m_fields(fields) {}

    void operator()(const amrex::Real time)
    {
        // cppcheck-suppress constVariableReference
        auto& field = m_fields.field;
        const auto& repo = field.repo();
        const int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& field_arrs = field(lev).arrays();
             const bool has_terrain = repo.int_field_exists("terrain_blank");
            if (has_terrain) {
	      amrex::ParmParse pp("DragTempForcing");
	      amrex::Real m_soil_temperature= 300.0;
	      pp.query("soil_temperature", m_soil_temperature);
                const IntField& m_terrain_blank =
                    repo.get_int_field("terrain_blank");
                const auto& blank_arrs = m_terrain_blank(lev).arrays();
                amrex::ParallelFor(
                    field(lev), [=] AMREX_GPU_DEVICE(
                                    int nbx, int i, int j, int k) noexcept {
                        field_arrs[nbx](i, j, k) =
                            (blank_arrs[nbx](i, j, k) == 1)
                                ? m_soil_temperature
                                : field_arrs[nbx](i, j, k);
                    });
            }
        }
        amrex::Gpu::streamSynchronize();
        m_fields.field.fillpatch(time);
    }

    PDEFields& m_fields;
};


} // namespace amr_wind::pde

#endif /* Temperature_OPS_H */
