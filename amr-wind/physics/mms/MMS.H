#ifndef MMS_H
#define MMS_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::mms {

using FuncDef = amrex::Real (*)(amrex::Real, amrex::Real, amrex::Real);

/** Method of manufactured solutions physics
 *  \ingroup physics
 *
 */
class MMS : public Physics::Register<MMS>
{
public:
    static std::string identifier() { return "MMS"; }

    explicit MMS(const CFDSim& sim);

    ~MMS() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void fill_src();

    amrex::Real compute_error(const int comp, const Field& field, FuncDef f);

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_density;
    Field& m_mms_vel_source;

    const int m_w = 18;
    const std::string m_output_fname = "mms.log";

    std::unordered_map<std::string, amrex::Real> m_params_map = {
        {"a", 1.05_rt},   {"b", 2.15_rt},    {"c", -3.2_rt},    {"d", 10.1_rt},
        {"beta", 2.2_rt}, {"gamma", 2.4_rt}, {"delta", 2.0_rt}, {"kx", 1.0_rt},
        {"ky", 1.0_rt},   {"kz", 1.0_rt},    {"nu", 1.0_rt}};
};
} // namespace amr_wind::mms

#endif /* MMS_H */
