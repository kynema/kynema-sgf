#ifndef FREE_SURFACE_DAMPING_H
#define FREE_SURFACE_DAMPING_H

#include "src/core/Physics.H"
#include "src/core/FieldRepo.H"
#include "src/core/SimTime.H"
#include "AMReX_Geometry.H"

namespace kynema_sgf {

class FreeSurfaceDamping : public Physics::Register<FreeSurfaceDamping>
{
public:
    static std::string identifier() { return "FreeSurfaceDamping"; }

    explicit FreeSurfaceDamping(CFDSim& /*sim*/);

    ~FreeSurfaceDamping() override;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override {
    };

    void post_init_actions() override {};

    void post_regrid_actions() override {};

    void pre_advance_work() override {};

    void post_advance_work() override;

private:
    FieldRepo& m_repo;

    bool m_vertical_velocity_damping{true};
    bool m_volume_fraction_damping{true};

    bool m_global_damping{true};

    amrex::Real m_time_scale{0.5_rt};

    amrex::Real m_length_xlo{0.0_rt};
    amrex::Real m_length_xhi{0.0_rt};
    amrex::Real m_length_ylo{0.0_rt};
    amrex::Real m_length_yhi{0.0_rt};

    amrex::Real m_rho1;
    amrex::Real m_rho2;
    amrex::Geometry m_geom;
    SimTime m_time;
};

} // namespace kynema_sgf

#endif /* FREE_SURFACE_DAMPING_H */
