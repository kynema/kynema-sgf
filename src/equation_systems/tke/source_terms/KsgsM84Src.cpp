#include <AMReX_Orientation.H>

#include "src/equation_systems/tke/source_terms/KsgsM84Src.H"
#include "src/CFDSim.H"
#include "src/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

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
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    const auto& repo = (this->m_tke).repo();
    const auto geom = repo.mesh().Geom(lev);

    const amrex::Real dx = geom.CellSize()[0];
    const amrex::Real dy = geom.CellSize()[1];
    const amrex::Real dz = geom.CellSize()[2];
    const amrex::Real ds = std::cbrt(dx * dy * dz);
    const amrex::Real Ceps = this->m_Ceps;
    const amrex::Real CepsGround = this->m_CepsGround;

    auto const& src_arrs = src_term.arrays();
    auto const& tlscale_arrs = m_turb_lscale(lev).const_arrays();
    auto const& shear_prod_arrs = m_shear_prod(lev).const_arrays();
    auto const& buoy_prod_arrs = m_buoy_prod(lev).const_arrays();
    auto const& dissip_arrs = m_dissip(lev).arrays();
    auto const& tke_arrs = m_tke(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            const auto& tlscale_arr = tlscale_arrs[nbx];
            const auto& shear_prod_arr = shear_prod_arrs[nbx];
            const auto& buoy_prod_arr = buoy_prod_arrs[nbx];
            const auto& dissip_arr = dissip_arrs[nbx];
            const auto& tke_arr = tke_arrs[nbx];

            dissip_arr(i, j, k) = calc_dissip(
                calc_ceps_local(Ceps, tlscale_arr(i, j, k), ds),
                tke_arr(i, j, k), tlscale_arr(i, j, k));
            src_arrs[nbx](i, j, k) += shear_prod_arr(i, j, k) +
                                      buoy_prod_arr(i, j, k) -
                                      dissip_arr(i, j, k);
        });

    // Wall boundary corrections via MFIter (boundary boxes are
    // MFIter-dependent)
    const auto& bctype = (this->m_tke).bc_type();
    for (amrex::MFIter mfi(src_term); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term.array(mfi);
        const auto& dissip_arr = m_dissip(lev).array(mfi);
        const auto& tke_arr = m_tke(lev).const_array(mfi);
        const auto& tlscale_arr = m_turb_lscale(lev).const_array(mfi);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            amrex::Orientation olo(dir, amrex::Orientation::low);
            if (bctype[olo] == BC::wall_model &&
                bx.smallEnd(dir) == geom.Domain().smallEnd(dir)) {
                amrex::Box blo = amrex::bdryLo(bx, dir, 1);
                if (blo.ok()) {
                    amrex::ParallelFor(
                        blo, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            src_arr(i, j, k) += dissip_arr(i, j, k);
                            dissip_arr(i, j, k) = calc_dissip(
                                CepsGround, tke_arr(i, j, k),
                                tlscale_arr(i, j, k));
                            src_arr(i, j, k) -= dissip_arr(i, j, k);
                        });
                } else {
                    amrex::Abort("Bad box extracted in KsgsM84Src");
                }
            }

            amrex::Orientation ohi(dir, amrex::Orientation::high);
            if (bctype[ohi] == BC::wall_model &&
                bx.bigEnd(dir) == geom.Domain().bigEnd(dir)) {
                amrex::Box bhi = amrex::bdryHi(bx, dir, 1);
                amrex::Abort(
                    "tke wall model is not supported on upper boundary");
                if (bhi.ok()) {
                    amrex::ParallelFor(
                        bhi, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            src_arr(i, j, k) += dissip_arr(i, j, k);
                            dissip_arr(i, j, k) = calc_dissip(
                                CepsGround, tke_arr(i, j, k),
                                tlscale_arr(i, j, k));
                            src_arr(i, j, k) -= dissip_arr(i, j, k);
                        });
                }
            }
        }
    }
}

} // namespace kynema_sgf::pde::tke
