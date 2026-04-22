#include <memory>

#include "amr-wind/physics/BoussinesqBubble.H"
#include "amr-wind/physics/BoussinesqBubbleFieldInit.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLWallFunction.H"

namespace amr_wind {

BoussinesqBubble::BoussinesqBubble(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register temperature equation
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");

    // Defer getting temperature field until PDE has been registered
    m_temperature = &(teqn.fields().field);

    // Instantiate the BoussinesqBubble field initializer
    m_field_init = std::make_unique<BoussinesqBubbleFieldInit>();
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::BoussinesqBubbleFieldInit
 */
void BoussinesqBubble::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& scalars = (*m_temperature)(level);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.tilebox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

} // namespace amr_wind
