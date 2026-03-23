#ifndef DAMBREAK_H
#define DAMBREAK_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/** Multiphase Dam break physics
 *  \ingroup multiphase_physics
 *
 */

namespace amr_wind {

class DamBreak : public Physics::Register<DamBreak>
{
public:
    static std::string identifier() { return "DamBreak"; }

    explicit DamBreak(CFDSim& sim);

    ~DamBreak() override = default;

    //! Initialize the levelset and velocity fields for Dam Break simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;
    Field& m_density;

    //! Initial DamBreak location
    amrex::Vector<amrex::Real> m_loc{0.0_rt, 0.0_rt, 0.0_rt};

    //! dam width value
    amrex::Real m_width{0.1_rt};

    //! dam height value
    amrex::Real m_height{0.25_rt};
};

} // namespace amr_wind

#endif /* DamBreak_H */
