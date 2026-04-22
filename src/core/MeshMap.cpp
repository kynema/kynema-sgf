#include "amr-wind/core/MeshMap.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

void MeshMap::declare_mapping_fields(const CFDSim& sim, int nghost)
{

    // declare nodal, cell-centered, and face-centered mesh mapping array
    m_mesh_scale_fac_cc = &(sim.repo().declare_cc_field(
        "mesh_scaling_factor_cc", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_nd = &(sim.repo().declare_nd_field(
        "mesh_scaling_factor_nd", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_xf = &(sim.repo().declare_xf_field(
        "mesh_scaling_factor_xf", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_yf = &(sim.repo().declare_yf_field(
        "mesh_scaling_factor_yf", AMREX_SPACEDIM, nghost, 1));
    m_mesh_scale_fac_zf = &(sim.repo().declare_zf_field(
        "mesh_scaling_factor_zf", AMREX_SPACEDIM, nghost, 1));

    // declare nodal, cell-centered, and face-centered mesh mapping detJ array
    m_mesh_scale_detJ_cc =
        &(sim.repo().declare_cc_field("mesh_scaling_detJ_cc", 1, nghost, 1));
    m_mesh_scale_detJ_nd =
        &(sim.repo().declare_nd_field("mesh_scaling_detJ_nd", 1, nghost, 1));
    m_mesh_scale_detJ_xf =
        &(sim.repo().declare_xf_field("mesh_scaling_detJ_xf", 1, nghost, 1));
    m_mesh_scale_detJ_yf =
        &(sim.repo().declare_yf_field("mesh_scaling_detJ_yf", 1, nghost, 1));
    m_mesh_scale_detJ_zf =
        &(sim.repo().declare_zf_field("mesh_scaling_detJ_zf", 1, nghost, 1));

    // declare nodal and cell-centered non-uniform mesh
    m_non_uniform_coord_cc = &(sim.repo().declare_cc_field(
        "non_uniform_coord_cc", AMREX_SPACEDIM, nghost, 1));
    m_non_uniform_coord_nd = &(sim.repo().declare_nd_field(
        "non_uniform_coord_nd", AMREX_SPACEDIM, nghost, 1));

    // TODO: Create BCNoOP fill patch operators for mesh scaling fields ?
}

} // namespace amr_wind
