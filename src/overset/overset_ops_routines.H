#ifndef OVERSET_OPS_ROUTINES_H
#define OVERSET_OPS_ROUTINES_H

#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFab.H"
#include "src/equation_systems/vof/volume_fractions.H"
#include "src/core/FieldRepo.H"

namespace kynema_sgf::overset_ops {

void iblank_to_mask(const IntField& iblank, IntField& maskf);

void iblank_node_to_mask_vof(
    const IntField& iblank, const Field& vof, IntField& maskf);

void prepare_mask_cell_for_mac(FieldRepo& repo);

void revert_mask_cell_after_mac(FieldRepo& repo);

// Swap pressure gradient values in overset region
void replace_gradp(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_gp0,
    const amrex::iMultiFab& mf_iblank);

// Apply pressure gradient to velocity field
void apply_pressure_gradient(
    amrex::MultiFab& mf_vel,
    const amrex::MultiFab& mf_density,
    const amrex::MultiFab& mf_gp,
    amrex::Real scaling_factor);

} // namespace kynema_sgf::overset_ops

#endif
