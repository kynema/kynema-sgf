#ifndef VORTEX_DIPOLE_H
#define VORTEX_DIPOLE_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
/** Vortex dipole Physics
 *   \ingroup physics
 */

class VortexDipole : public Physics::Register<VortexDipole>
{
public:
    static std::string identifier() { return "VortexDipole"; }

    explicit VortexDipole(const CFDSim& sim);

    ~VortexDipole() override = default;

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

    //! Initial location of the left vortex
    amrex::Vector<amrex::Real> m_loc_left{-0.1_rt, 0., 0.};

    //! Initial location of the right vortex
    amrex::Vector<amrex::Real> m_loc_right{0.1_rt, 0., 0.};

    //! Background velocity
    amrex::Vector<amrex::Real> m_bvel{0., 0., 0.};

    //! Initial vorticity
    amrex::Real m_omegaEmag{299.5283853752_rt};

    //! Vortex
    amrex::Real m_r0{0.1_rt};
};

} // namespace amr_wind

#endif /* VORTEX_DIPOLE_H */
