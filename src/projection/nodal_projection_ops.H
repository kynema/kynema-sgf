#ifndef NODAL_PROJECTION_OPS_H
#define NODAL_PROJECTION_OPS_H

#include "AMReX_BC_TYPES.H"
#include "src/core/field_ops.H"
#include "src/boundary_conditions/field_boundary_fill/FieldBoundary.H"
#include "src/wind_energy/ABL.H"
#include "src/ocean_waves/OceanWaves.H"

namespace kynema_sgf::nodal_projection {

void set_inflow_velocity(
    kynema_sgf::FieldBoundaryMgr::TypeVector& fbs,
    kynema_sgf::Field& vel_fld,
    int lev,
    amrex::Real time,
    amrex::MultiFab& vel_mfab,
    int nghost);

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    amrex::Orientation::Side side,
    kynema_sgf::Field& pressure,
    const amrex::Array<int, AMREX_SPACEDIM>& is_periodic);

void apply_dirichlet_vel(
    amrex::MultiFab& mf_velocity, amrex::iMultiFab& mf_iblank);

void enforce_inout_solvability(
    kynema_sgf::Field& velocity,
    const amrex::Vector<amrex::Geometry>& geom,
    int num_levels);

} // namespace kynema_sgf::nodal_projection

#endif
