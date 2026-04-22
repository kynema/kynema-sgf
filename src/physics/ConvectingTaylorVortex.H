#ifndef CONVECTINGTAYLORVORTEX_H
#define CONVECTINGTAYLORVORTEX_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ctv {

namespace {
struct UExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*u0*/,
        amrex::Real /*v0*/,
        amrex::Real /*omega*/,
        amrex::Real /*x*/,
        amrex::Real /*y*/,
        amrex::Real /*t*/) const;
    const int m_comp{0};
};

struct VExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*u0*/,
        amrex::Real /*v0*/,
        amrex::Real /*omega*/,
        amrex::Real /*x*/,
        amrex::Real /*y*/,
        amrex::Real /*t*/) const;
    const int m_comp{1};
};

struct WExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/) const;
    const int m_comp{2};
};

struct GpxExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*u0*/,
        amrex::Real /*unused*/,
        amrex::Real /*omega*/,
        amrex::Real /*x*/,
        amrex::Real /*unused*/,
        amrex::Real /*t*/) const;
    const int m_comp{0};
};

struct GpyExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*unused*/,
        amrex::Real /*v0*/,
        amrex::Real /*omega*/,
        amrex::Real /*unused*/,
        amrex::Real /*y*/,
        amrex::Real /*t*/) const;
    const int m_comp{1};
};

struct GpzExact
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/) const;
    const int m_comp{2};
};

} // namespace

/** Convecting Taylor Vortex physics
 *  \ingroup physics
 */
class ConvectingTaylorVortex : public Physics::Register<ConvectingTaylorVortex>
{
public:
    static std::string identifier() { return "ConvectingTaylorVortex"; }

    explicit ConvectingTaylorVortex(const CFDSim& sim);

    ~ConvectingTaylorVortex() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    template <typename T>
    amrex::Real compute_error(const Field& /*field*/);

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

    [[nodiscard]] amrex::Real get_u0() const { return m_u0; }

    [[nodiscard]] amrex::Real get_v0() const { return m_v0; }

    [[nodiscard]] amrex::Real get_omega() const { return m_omega; }

private:
    const amr_wind::SimTime& m_time;
    const amr_wind::CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    Field& m_velocity;
    Field& m_gradp;
    Field& m_density;

    void output_error();

    //! initial density value
    amrex::Real m_rho{1.0_rt};

    //! mean x-velocity value
    amrex::Real m_u0{1.0_rt};

    //! mean y-velocity value
    amrex::Real m_v0{1.0_rt};

    //! damping
    amrex::Real m_omega;

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"ctv.log"};

    bool m_activate_pressure{false};
    bool m_mesh_mapping{false};
};
} // namespace amr_wind::ctv

#endif /* ConvectingTaylorVortex_H */
