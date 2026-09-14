#ifndef FLATHER_H
#define FLATHER_H

#include "src/core/Field.H"
#include "src/core/IntField.H"
#include "src/CFDSim.H"
#include "src/boundary_conditions/field_boundary_fill/FieldBoundary.H"
#include "src/utilities/MultiLevelVector.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_REAL.H"
#include <limits>

using namespace amrex::literals;

namespace kynema_sgf {

/** Flather-style open boundary for velocity
 *
 *  This boundary reads current boundary ghost-cell values and adjacent interior
 *  values, computes a blended target, and writes updated boundary values.
 */
class Flather : public FieldBoundary::Register<Flather>
{
public:
    static std::string identifier() { return "Flather"; }

    explicit Flather(CFDSim& sim);

    void post_init_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override {}

    void set_velocity(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab,
        int dcomp = 0,
        int orig_comp = 0) const override;

    void update_flather_variables(
        const int lev,
        const FieldState fstate,
        const bool use_mac_fields = false)
    {
        compute_boundary_z_averages(lev, fstate, use_mac_fields);
    }

    void accumulate_boundary(
        int current_level,
        int idir,
        int phase_switch,
        bool is_low,
        MultiLevelVector& out_uvec,
        MultiLevelVector& out_hvec,
        bool sample_boundary,
        FieldState fstate,
        bool use_mac_fields = false) const;

private:
    void compute_internal_z_averages();
    void compute_boundary_z_averages(
        int lev, FieldState fstate, bool use_mac_fields = false);

    const CFDSim& m_sim;
    const kynema_sgf::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    const Field& m_u_mac;
    const Field& m_v_mac;
    const Field& m_vof;
    const IntField* m_terrain_blank{nullptr};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
    amrex::Real m_velocity_scale_factor_max{2.0_rt};
    amrex::Real m_velocity_scale_factor_min{0.25_rt};

    MultiLevelVector m_xlo_uhliq{FieldLoc::CELL};
    MultiLevelVector m_xhi_uhliq{FieldLoc::CELL};
    MultiLevelVector m_ylo_uhliq{FieldLoc::CELL};
    MultiLevelVector m_yhi_uhliq{FieldLoc::CELL};
    MultiLevelVector m_xlo_uhmix{FieldLoc::CELL};
    MultiLevelVector m_xhi_uhmix{FieldLoc::CELL};
    MultiLevelVector m_ylo_uhmix{FieldLoc::CELL};
    MultiLevelVector m_yhi_uhmix{FieldLoc::CELL};
    MultiLevelVector m_xlo_bnd_uh{FieldLoc::CELL};
    MultiLevelVector m_xhi_bnd_uh{FieldLoc::CELL};
    MultiLevelVector m_ylo_bnd_uh{FieldLoc::CELL};
    MultiLevelVector m_yhi_bnd_uh{FieldLoc::CELL};
    MultiLevelVector m_xlo_h{FieldLoc::CELL};
    MultiLevelVector m_xhi_h{FieldLoc::CELL};
    MultiLevelVector m_ylo_h{FieldLoc::CELL};
    MultiLevelVector m_yhi_h{FieldLoc::CELL};
    MultiLevelVector m_xlo_bnd_h{FieldLoc::CELL};
    MultiLevelVector m_xhi_bnd_h{FieldLoc::CELL};
    MultiLevelVector m_ylo_bnd_h{FieldLoc::CELL};
    MultiLevelVector m_yhi_bnd_h{FieldLoc::CELL};
};

} // namespace kynema_sgf

#endif /* FLATHER_H */
