#include "amr-wind/equation_systems/icns/source_terms/SynthTurbForcing.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

/** Synthetic Turbulence forcing term
 *
 *
 *
 */
SynthTurbForcing::SynthTurbForcing(const CFDSim& sim)
    : m_turb_force(sim.repo().get_field("synth_turb_forcing"))
{
    if (!sim.physics_manager().contains("SyntheticTurbulence")) {
        amrex::Abort(
            "SynthTurbForcing: SyntheticTurbulence physics not enabled");
    }
}

SynthTurbForcing::~SynthTurbForcing() = default;

void SynthTurbForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& turb_force_arr = m_turb_force(lev).array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        src_term(i, j, k, 0) += turb_force_arr(i, j, k, 0);
        src_term(i, j, k, 1) += turb_force_arr(i, j, k, 1);
        src_term(i, j, k, 2) += turb_force_arr(i, j, k, 2);
    });
}

} // namespace amr_wind::pde::icns
