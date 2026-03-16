#ifndef TAYLORGREENVORTEX_H
#define TAYLORGREENVORTEX_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/** Taylor Green Vortex physics
 *  \ingroup physics
 */
class TaylorGreenVortex : public Physics::Register<TaylorGreenVortex>
{
public:
    static std::string identifier() { return "TaylorGreenVortex"; }

    explicit TaylorGreenVortex(const CFDSim& sim);

    ~TaylorGreenVortex() override = default;
    //! Initialize the temperature and velocity fields for BoussinesqBubble
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;

    //! initial density value
    amrex::Real m_rho{1.0_rt};
};

} // namespace amr_wind

#endif /* TaylorGreenVortex_H */
