#ifndef NODAL_PROJECTION_OPS_H
#define NODAL_PROJECTION_OPS_H

#include "AMReX_BC_TYPES.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/ocean_waves/OceanWaves.H"

namespace amr_wind::nodal_projection {

void set_inflow_velocity(
    amr_wind::PhysicsMgr& phy_mgr,
    amr_wind::Field& vel_fld,
    int lev,
    amrex::Real time,
    amrex::MultiFab& vel_mfab,
    int nghost);

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    amrex::Orientation::Side side,
    amr_wind::Field& pressure,
    const amrex::Array<int, AMREX_SPACEDIM>& is_periodic);

void apply_dirichlet_vel(
    amrex::MultiFab& mf_velocity, amrex::iMultiFab& mf_iblank);

void enforce_inout_solvability(
    amr_wind::Field& velocity,
    const amrex::Vector<amrex::Geometry>& geom,
    int num_levels);

} // namespace amr_wind::nodal_projection

#endif
