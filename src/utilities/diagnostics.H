#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include "AMReX_MultiFabUtil.H"

namespace amr_wind::diagnostics {

void make_mask_addend(
    amrex::MultiFab& mfab,
    const amrex::MultiFab& mfab_mask,
    amrex::Real mask_val,
    amrex::Real set_val);

void get_field_extrema(
    amrex::Real& field_max_val,
    amrex::Real& field_min_val,
    const amr_wind::Field& field,
    int comp,
    int ncomp,
    int nghost);

bool get_field_extrema(
    amrex::Real& field_max_val,
    amrex::Real& field_min_val,
    const amr_wind::Field& field,
    const amr_wind::Field& field_mask,
    amrex::Real mask_val,
    int comp,
    int ncomp,
    int nghost);

amrex::Real get_vel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    int vdir,
    amrex::Real factor);

amrex::Real get_vel_max(
    const amrex::MultiFab& vel, const amrex::iMultiFab& level_mask, int vdir);

amrex::Real get_vel_min(
    const amrex::MultiFab& vel, const amrex::iMultiFab& level_mask, int vdir);

amrex::Real get_vel_loc(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    int vdir,
    int ldir,
    amrex::Real vel_max,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx);

amrex::Real get_macvel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    int vdir,
    amrex::Real factor);

amrex::Real get_macvel_max(
    const amrex::MultiFab& vel, const amrex::iMultiFab& level_mask, int vdir);

amrex::Real get_macvel_min(
    const amrex::MultiFab& vel, const amrex::iMultiFab& level_mask, int vdir);

amrex::Real get_macvel_loc(
    const amrex::MultiFab& macvel,
    const amrex::iMultiFab& level_mask,
    int vdir,
    int ldir,
    amrex::Real vel_max,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx);

amrex::Array<amrex::Real, 24> PrintMaxVelLocations(
    const amr_wind::FieldRepo& repo, const std::string& header);

amrex::Array<amrex::Real, 24> PrintMaxMACVelLocations(
    const amr_wind::FieldRepo& repo, const std::string& header);

} // namespace amr_wind::diagnostics

#endif
