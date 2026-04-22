#ifndef EKMANSPIRAL_H
#define EKMANSPIRAL_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/** Ekman Spiral physics
 *  \ingroup physics
 */
class EkmanSpiral : public Physics::Register<EkmanSpiral>
{
public:
    static std::string identifier() { return "EkmanSpiral"; }

    explicit EkmanSpiral(const CFDSim& sim);

    ~EkmanSpiral() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    template <typename T>
    amrex::Real compute_error(const Field& /*field*/);

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_density;

    void output_error();

    //! initial density value
    amrex::Real m_rho{1.0_rt};
    //! Ekman depth
    amrex::Real m_DE{0.0_rt};
    //! velocity magnitude x-direction only
    amrex::Real m_vel{15.0_rt};

    //! output precision
    const int m_w{18};

    //! error log file
    const std::string m_output_fname{"ekman_spiral.log"};
};
} // namespace amr_wind

#endif /* EkmanSpiral_H */
