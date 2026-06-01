#ifndef OCEANWAVESBOUNDARY_H
#define OCEANWAVESBOUNDARY_H

#include "src/core/Field.H"
#include "src/core/IntField.H"
#include "src/CFDSim.H"
#include "src/boundary_conditions/field_boundary_fill/FieldBoundary.H"
#include "AMReX_Gpu.H"
#include <AMReX_BndryRegister.H>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

/** Interface for Ocean Waves boundary data
 *
 *  This class performs the necessary field operations to provide
 *  available wave data at the boundaries
 *
 */
class OceanWavesBoundary : public FieldBoundary::Register<OceanWavesBoundary>
{
public:
    static std::string identifier() { return "OceanWavesBoundary"; }

    explicit OceanWavesBoundary(CFDSim& /*sim*/);

    void post_init_actions() override;

    void pre_advance_work() override;

    void pre_predictor_work() override;

    void post_advance_work() override {}

    void set_velocity(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab,
        int dcomp = 0,
        int orig_comp = 0) const override;

    void
    set_vof(int lev, amrex::Real time, const Field& fld, amrex::MultiFab& mfab)
        const;

    void set_density(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

    void set_inflow_sibling_velocity(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs) const;

    void record_boundary_data_time(const amrex::Real time)
    {
        m_bndry_time = time;
    }

    [[nodiscard]] amrex::Real get_boundary_data_time() const
    {
        return m_bndry_time;
    }

private:
    const kynema_sgf::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field* m_ow_velocity{nullptr};
    Field* m_ow_vof{nullptr};
    IntField* m_terrain_blank_ptr{nullptr};

    // Check for single-phase simulation
    bool m_vof_exists{true};

    // Check for terrain
    bool m_terrain_exists{false};

    // Store time corresponding to current boundary data
    amrex::Real m_bndry_time{0.0_rt};

    // Will only remain negative if not needed
    amrex::Real m_rho1{-1.0_rt};
};

} // namespace kynema_sgf

#endif /* OCEANWAVESBOUNDARY_H */
