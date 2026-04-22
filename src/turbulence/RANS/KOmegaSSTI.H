#ifndef KOMEGASSTI_H
#define KOMEGASSTI_H

#include "amr-wind/turbulence/RANS/KOmegaSST.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/sdr/SDR.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
KOmegaSST<Transport>::KOmegaSST(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_f1(sim.repo().declare_field("f1", 1, 1, 1))
    , m_shear_prod(sim.repo().declare_field("shear_prod", 1, 1, 1))
    , m_buoy_term(sim.repo().declare_field("buoyancy_term", 1, 1, 1))
    , m_diss(sim.repo().declare_field("dissipation", 1, 1, 1))
    , m_sdr_src(sim.repo().declare_field("omega_src", 1, 1, 1))
    , m_sdr_diss(sim.repo().declare_field("sdr_dissipation", 1, 1, 1))
    , m_rho(sim.repo().get_field("density"))
    , m_walldist(sim.repo().declare_field("wall_dist", 1, 1, 1))
{
    m_walldist.set_default_fillpatch_bc(sim.time());
    m_walldist.fillpatch_on_regrid() = true;

    auto& tke_eqn =
        sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);

    auto& sdr_eqn =
        sim.pde_manager().register_transport_pde(pde::SDR::pde_name());
    m_sdr = &(sdr_eqn.fields().field);

    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(pde::TKE::pde_name(), {"KwSSTSrc"});
}

template <typename Transport>
KOmegaSST<Transport>::~KOmegaSST() = default;

} // namespace amr_wind::turbulence

#endif /* KOMEGASSTI_H */
