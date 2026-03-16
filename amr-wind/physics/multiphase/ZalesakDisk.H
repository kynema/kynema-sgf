#ifndef ZALESAKDISK_H
#define ZALESAKDISK_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/** Multiphase Zalesak disk physics
 *  \ingroup multiphase_physics
 *
 */

namespace amr_wind {

class ZalesakDisk : public Physics::Register<ZalesakDisk>
{

    static_assert(
        AMREX_SPACEDIM == 3, "ZalesakDisk requires 3 dimensional mesh");

public:
    static std::string identifier() { return "ZalesakDisk"; }

    explicit ZalesakDisk(CFDSim& sim);

    ~ZalesakDisk() override = default;

    //! Initialize the levelset and velocity fields for Vortex Patch simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override;

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;
    Field& m_density;

    //! Initial ZalesakDisk location
    amrex::Vector<amrex::Real> m_loc{0.5_rt, 0.72_rt, 0.24_rt};

    //! sphere radius value
    amrex::Real m_radius{0.16_rt};

    //! slot half width value
    amrex::Real m_halfwidth{0.04_rt};

    //! slot depth
    amrex::Real m_depth{0.2_rt};

    //! prescribed velocity period
    amrex::Real m_TT{6.0_rt};
};

} // namespace amr_wind

#endif /* ZalesakDisk_H */
