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

//! Multiply the face diffusion coefficients by the no-slip wall-flux factors
//! declared by ImmersedTerrain. No-op when the factor fields do not exist.
void apply_immersed_interface(
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fc,
    const kynema_sgf::FieldRepo& repo,
    int lev);

//! Effective density rho (1 + beta C dt) for the implicit diffusion solve
//! when the immersed terrain drag is applied implicitly (ImmersedTerrain
//! declares terrain_drag_rate). Pins the body cells during the solve so that
//! they are not filled from the fluid and then emptied by the projection.
//! Returns nullptr when the field does not exist.
std::unique_ptr<kynema_sgf::ScratchField> immersed_effective_density(
    const kynema_sgf::FieldRepo& repo,
    const kynema_sgf::Field& density,
    amrex::Real dt);

void viscosity_to_uniform_space(
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& b,
    const kynema_sgf::FieldRepo& repo,
    int lev);

} // namespace diffusion

#endif /* DIFFUSION_H */
