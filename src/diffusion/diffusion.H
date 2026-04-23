#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "AMReX_MLMG.H"
#include "src/core/FieldRepo.H"

namespace diffusion {

amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
get_diffuse_tensor_bc(
    kynema_sgf::Field& velocity, amrex::Orientation::Side side);

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>
get_diffuse_scalar_bc(kynema_sgf::Field& scalar, amrex::Orientation::Side side);

amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> average_velocity_eta_to_faces(
    const amrex::Geometry& geom, amrex::MultiFab const& cc_eta);

void fixup_eta_on_domain_faces(
    const amrex::Geometry& geom,
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fc,
    amrex::MultiFab const& cc);

void viscosity_to_uniform_space(
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& b,
    const kynema_sgf::FieldRepo& repo,
    int lev);

} // namespace diffusion

#endif /* DIFFUSION_H */
