#include <memory>

#include "src/physics/RayleighTaylor.H"
#include "src/physics/RayleighTaylorFieldInit.H"
#include "src/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

RayleighTaylor::RayleighTaylor(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_field_init(std::make_unique<RayleighTaylorFieldInit>())
{}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 *
 *  \sa kynema_sgf::RayleighTaylorFieldInit
 */
void RayleighTaylor::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    velocity.setVal(0.0_rt);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.tilebox();

        (*m_field_init)(vbx, geom, density.array(mfi));
    }
}

} // namespace kynema_sgf
