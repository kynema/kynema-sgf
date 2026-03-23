#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KsgsM84Src.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::tke {

KsgsM84Src::KsgsM84Src(const CFDSim& sim)
    : m_turb_lscale(sim.repo().get_field("turb_lscale"))
    , m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_buoy_prod(sim.repo().get_field("buoy_prod"))
    , m_dissip(sim.repo().get_field("dissipation"))
    , m_tke(sim.repo().get_field("tke"))
{
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "OneEqKsgsM84");
    auto coeffs = sim.turbulence_model().model_coeffs();
    m_Ceps = coeffs["Ceps"];
    m_CepsGround = (3.9_rt / 0.93_rt) * m_Ceps;
}

KsgsM84Src::~KsgsM84Src() = default;

void KsgsM84Src::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
    const auto& dissip_arr = (this->m_dissip)(lev).array(mfi);
    const auto& tke_arr = (this->m_tke)(lev).array(mfi);
    const amrex::Real Ceps = this->m_Ceps;
    const amrex::Real CepsGround = this->m_CepsGround;

    const auto& repo = (this->m_tke).repo();
    const auto geom = repo.mesh().Geom(lev);

    const amrex::Real dx = geom.CellSize()[0];
    const amrex::Real dy = geom.CellSize()[1];
    const amrex::Real dz = geom.CellSize()[2];
    const amrex::Real ds = std::cbrt(dx * dy * dz);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        dissip_arr(i, j, k) = calc_dissip(
            calc_ceps_local(Ceps, tlscale_arr(i, j, k), ds), tke_arr(i, j, k),
            tlscale_arr(i, j, k));
        src_term(i, j, k) += shear_prod_arr(i, j, k) + buoy_prod_arr(i, j, k) -
                             dissip_arr(i, j, k);
    });

    const auto& bctype = (this->m_tke).bc_type();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::Orientation olo(dir, amrex::Orientation::low);
        if (bctype[olo] == BC::wall_model &&
            bx.smallEnd(dir) == geom.Domain().smallEnd(dir)) {
            amrex::Box blo = amrex::bdryLo(bx, dir, 1);
            if (blo.ok()) {
                amrex::ParallelFor(
                    blo, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        // Remove previous dissipation term
                        src_term(i, j, k) += dissip_arr(i, j, k);
                        // Calculate new dissipation term
                        dissip_arr(i, j, k) = calc_dissip(
                            CepsGround, tke_arr(i, j, k), tlscale_arr(i, j, k));
                        // Add new dissiaption term
                        src_term(i, j, k) -= dissip_arr(i, j, k);
                    });
            } else {
                amrex::Abort("Bad box extracted in KsgsM84Src");
            }
        }

        amrex::Orientation ohi(dir, amrex::Orientation::high);
        if (bctype[ohi] == BC::wall_model &&
            bx.bigEnd(dir) == geom.Domain().bigEnd(dir)) {
            amrex::Box bhi = amrex::bdryHi(bx, dir, 1);
            amrex::Abort("tke wall model is not supported on upper boundary");
            if (bhi.ok()) {
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        // Remove previous dissipation term
                        src_term(i, j, k) += dissip_arr(i, j, k);
                        // Calculate new dissipation term
                        dissip_arr(i, j, k) = calc_dissip(
                            CepsGround, tke_arr(i, j, k), tlscale_arr(i, j, k));
                        // Add new dissiaption term
                        src_term(i, j, k) -= dissip_arr(i, j, k);
                    });
            }
        }
    }
}

} // namespace amr_wind::pde::tke
