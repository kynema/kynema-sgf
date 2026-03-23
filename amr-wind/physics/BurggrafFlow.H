#ifndef BURGGRAFFLOW_H
#define BURGGRAFFLOW_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::burggraf {

namespace {

struct UExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real /*x*/, amrex::Real /*y*/) const;
    const int m_comp{0};
};

struct VExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real /*x*/, amrex::Real /*y*/) const;
    const int m_comp{1};
};

struct WExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real /*x*/, amrex::Real /*y*/) const;
    const int m_comp{2};
};

} // namespace

/** Burggraf Flow physics
 *  \defgroup burggraf_flow Burggraf flow
 *  \ingroup physics
 */

class BurggrafFlow : public Physics::Register<BurggrafFlow>
{
public:
    static std::string identifier() { return "BurggrafFlow"; }

    explicit BurggrafFlow(const CFDSim& sim);

    ~BurggrafFlow() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    template <typename T>
    amrex::Real compute_error(const Field& /*field*/);

    void apply_lid_boundary_condition();

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {};

    void post_advance_work() override;

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    Field& m_velocity;
    Field& m_density;
    Field& m_source;

    void output_error();

    //! initial density value
    amrex::Real m_rho{1.0_rt};

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"burggraf.log"};
};
} // namespace amr_wind::burggraf

#endif
