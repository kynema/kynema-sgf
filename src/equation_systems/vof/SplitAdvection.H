#ifndef SPLITADVECTION_H
#define SPLITADVECTION_H

#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::multiphase {

void split_advection_step(
    int isweep,
    int iorder,
    int nlevels,
    Field& dof_field,
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> const& fluxes,
    ScratchField& fluxC,
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> const& advas,
    Field const& u_mac,
    Field const& v_mac,
    Field const& w_mac,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> BCs,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real time,
    amrex::Real dt,
    bool rm_debris);

void split_compute_fluxes(
    int lev,
    amrex::Box const& bx,
    int isweep,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::Array4<amrex::Real> const& aax,
    amrex::Array4<amrex::Real> const& aay,
    amrex::Array4<amrex::Real> const& aaz,
    amrex::Array4<amrex::Real> const& fx,
    amrex::Array4<amrex::Real> const& fy,
    amrex::Array4<amrex::Real> const& fz,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> BCs,
    amrex::Real* p,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt);

void split_compute_sum(
    int lev,
    amrex::Box const& bx,
    int isweep,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxC,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt);

void cmask_loop(
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxC);

void debris_loop(
    amrex::Box const& bx, amrex::Array4<amrex::Real> const& volfrac);

void sweep_fluxes(
    int dir,
    amrex::Box const& bx,
    amrex::Real dtdx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real> const& vofL,
    amrex::Array4<amrex::Real> const& vofR);

void sweep_balance(
    int dir,
    amrex::Real dt,
    amrex::Real dxi,
    amrex::Box const& bx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxF,
    amrex::Array4<amrex::Real const> const& fluxC);

} // namespace amr_wind::multiphase
#endif /* SplitAdvection_H */
