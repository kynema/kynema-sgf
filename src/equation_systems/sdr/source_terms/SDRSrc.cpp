#include <AMReX_Orientation.H>

#include "src/equation_systems/sdr/source_terms/SDRSrc.H"
#include "src/CFDSim.H"
#include "src/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

SDRSrc::SDRSrc(const CFDSim& sim)
    : m_sdr_src(sim.repo().get_field("omega_src"))
    , m_sdr_diss(sim.repo().get_field("sdr_dissipation"))
{}

SDRSrc::~SDRSrc() = default;

void SDRSrc::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& sdr_src_arrs = (this->m_sdr_src)(lev).const_arrays();
    const auto& sdr_diss_arrs = (this->m_sdr_diss)(lev).const_arrays();

    const amrex::Real factor = (fstate == FieldState::NPH) ? 0.5_rt : 1.0_rt;

    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            src_arrs[nbx](i, j, k) += (factor * sdr_diss_arrs[nbx](i, j, k)) +
                                      sdr_src_arrs[nbx](i, j, k);
        });
}

} // namespace kynema_sgf::pde::tke
