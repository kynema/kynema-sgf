#include "amr-wind/physics/HybridRANSLESABL.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind::hybrid_rans_les_abl {

HybridRANSLESABL::HybridRANSLESABL(const CFDSim& sim) : m_sim(sim) {}

/** Initialize the sdr field at the beginning of the simulation.
 */
void HybridRANSLESABL::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    const amrex::Real Ce = this->m_Ce;
    const amrex::Real dx = geom.CellSize()[0];
    const amrex::Real dy = geom.CellSize()[1];
    const amrex::Real dz = geom.CellSize()[2];
    const amrex::Real ds = std::cbrt(dx * dy * dz);

    const auto& tke_arrs = (*m_tke)(level).const_arrays();
    const auto& sdr_arrs = (*m_sdr)(level).arrays();

    amrex::ParallelFor(
        (*m_tke)(level), m_tke->num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            sdr_arrs[nbx](i, j, k) =
                std::sqrt(tke_arrs[nbx](i, j, k)) / (Ce * ds);
        });
    amrex::Gpu::streamSynchronize();
}

void HybridRANSLESABL::pre_init_actions()
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_init_actions");

    AMREX_ALWAYS_ASSERT(
        m_sim.turbulence_model().model_name() == "OneEqKsgsM84");
    auto t_coeffs = m_sim.turbulence_model().model_coeffs();
    m_Ce = t_coeffs["Ce"];

    m_tke = &(m_sim.repo().get_field("tke"));

    if (m_sim.repo().field_exists("sdr")) {
        m_sdr = &(m_sim.repo().get_field("sdr"));
    } else {
        m_sdr =
            &(m_sim.repo().declare_field("sdr", 1, m_tke->num_grow()[0], 1));
    }
}

void HybridRANSLESABL::post_init_actions()
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");

    compute_sdr_impl();
}

void HybridRANSLESABL::post_advance_work()
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::post_advance_work");

    compute_sdr_impl();
}

// Update sdr field based on sfs ke
void HybridRANSLESABL::compute_sdr_impl()
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::compute_sdr_impl");

    auto* tke = this->m_tke;
    auto* sdr = this->m_sdr;
    const amrex::Real Ce = this->m_Ce;

    const auto& geom_vec = m_sim.repo().mesh().Geom();
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        const auto& tke_arrs = (*tke)(lev).const_arrays();
        const auto& sdr_arrs = (*sdr)(lev).arrays();

        amrex::ParallelFor(
            (*tke)(lev), tke->num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                sdr_arrs[nbx](i, j, k) =
                    std::sqrt(tke_arrs[nbx](i, j, k)) / (Ce * ds);
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind::hybrid_rans_les_abl
