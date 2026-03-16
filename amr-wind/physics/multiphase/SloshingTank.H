#ifndef SLOSHINGTANK_H
#define SLOSHINGTANK_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/** Multiphase Sloshing Tank 3D physics
 *  \ingroup multiphase_physics
 *
 */

namespace amr_wind {

class SloshingTank : public Physics::Register<SloshingTank>
{
public:
    static std::string identifier() { return "SloshingTank"; }

    explicit SloshingTank(CFDSim& sim);

    ~SloshingTank() override = default;

    //! Initialize the levelset and velocity fields for Sloshing Tank
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_levelset;
    Field& m_pressure;

    //! Initial free surface amplitude magnitude
    amrex::Real m_amplitude{0.1_rt};

    //! Peak enhancement factor
    amrex::Real m_kappa{0.25_rt};

    //! Initial zero-level free-surface water depth
    amrex::Real m_waterlevel{0.0_rt};

    bool m_init_p{false};

    //! Stuff to get from MultiPhase physics
    amrex::Real m_rho1{1000.};
    amrex::Real m_rho2{1.};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
};

} // namespace amr_wind
#endif /* SloshingTank_H */
