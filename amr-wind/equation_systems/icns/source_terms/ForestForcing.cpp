#include "amr-wind/equation_systems/icns/source_terms/ForestForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind::pde::icns {

ForestForcing::ForestForcing(const CFDSim& sim)
    : m_sim(sim), m_velocity(sim.repo().get_field("velocity"))
{}

ForestForcing::~ForestForcing() = default;

void ForestForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const bool has_forest = this->m_sim.repo().field_exists("forest_drag");
    if (!has_forest) {
        amrex::Abort("Need a forest to use this source term");
    }
    auto* const m_forest_drag = &this->m_sim.repo().get_field("forest_drag");
    const auto& forest_drag = (*m_forest_drag)(lev).const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real ux = vel(i, j, k, 0);
        const amrex::Real uy = vel(i, j, k, 1);
        const amrex::Real uz = vel(i, j, k, 2);
        const amrex::Real windspeed =
            std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
        src_term(i, j, k, 0) -= forest_drag(i, j, k) * ux * windspeed;
        src_term(i, j, k, 1) -= forest_drag(i, j, k) * uy * windspeed;
        src_term(i, j, k, 2) -= forest_drag(i, j, k) * uz * windspeed;
    });
}

} // namespace amr_wind::pde::icns
