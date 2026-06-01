#include "src/equation_systems/icns/source_terms/ForestForcing.H"
#include "src/utilities/IOManager.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "src/wind_energy/ABL.H"

namespace kynema_sgf::pde::icns {

ForestForcing::ForestForcing(const CFDSim& sim)
    : m_sim(sim), m_velocity(sim.repo().get_field("velocity"))
{}

ForestForcing::~ForestForcing() = default;

void ForestForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const bool has_forest = this->m_sim.repo().field_exists("forest_drag");
    if (!has_forest) {
        amrex::Abort("Need a forest to use this source term");
    }
    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    auto const& forest_arrs =
        this->m_sim.repo().get_field("forest_drag")(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            const amrex::Real ux = vel_arrs[nbx](i, j, k, 0);
            const amrex::Real uy = vel_arrs[nbx](i, j, k, 1);
            const amrex::Real uz = vel_arrs[nbx](i, j, k, 2);
            const amrex::Real windspeed =
                std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
            src_arrs[nbx](i, j, k, n) -= forest_arrs[nbx](i, j, k) *
                                         vel_arrs[nbx](i, j, k, n) * windspeed;
        });
}

} // namespace kynema_sgf::pde::icns
