#ifndef SCALARADVECTION_H
#define SCALARADVECTION_H

#include <numbers>
#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/** Scalar Advection physics
 *  \ingroup physics
 *
 *  \sa ScalarAdvectionFieldInit
 */
class ScalarAdvection : public Physics::Register<ScalarAdvection>
{
public:
    static std::string identifier() { return "ScalarAdvection"; }

    explicit ScalarAdvection(CFDSim& sim);

    ~ScalarAdvection() override = default;

    //! Initialize the fields for ScalarAdvection simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    template <typename Shape>
    amrex::Vector<amrex::Real> compute_error(const Shape& scalar_function);

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

    template <typename Shape>
    void initialize_scalar(const Shape& scalar_function);

private:
    const SimTime& m_time;
    const FieldRepo& m_repo;
    Field& m_velocity;
    Field& m_density;
    std::string m_output_fname = "error.log";
    const int m_w = 20;

    //! Scalar field
    Field* m_scalar{nullptr};

    //! Flow velocity in x-direction
    amrex::Real m_u{1.0_rt};
    amrex::Real m_v{0.0_rt};

    //! Initial wave/pulse packet location
    amrex::Real m_x0{0.1_rt};
    amrex::Real m_y0{0.1_rt};

    //! Initial wave/pulse amplitude
    amrex::Real m_amplitude{1.0_rt};

    //! Initial wave/pulse width (standard deviation for GaussianPulse and
    //! GaussianWavePacket)
    amrex::Real m_x_width{0.01_rt};
    amrex::Real m_y_width{0.01_rt};

    //! Wave numbers
    amrex::Real m_x_wavenumber{16.0_rt * std::numbers::pi_v<amrex::Real>};
    amrex::Real m_y_wavenumber{16.0_rt * std::numbers::pi_v<amrex::Real>};

    //! Initial wave/pulse shape
    std::string m_shape{"gaussian"};

    //! Density
    amrex::Real m_rho{1.0_rt};
};

struct GaussianWavePacket
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*x*/,
        amrex::Real /*x0*/,
        amrex::Real /*amplitude*/,
        amrex::Real /*x_width*/,
        amrex::Real /*x_wavenumber*/) const;
};

struct GaussianPulseFV
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*x*/,
        amrex::Real /*unused*/,
        amrex::Real /*dx*/,
        amrex::Real /*unused*/,
        amrex::Real /*x0*/,
        amrex::Real /*unused*/,
        amrex::Real /*amplitude*/,
        amrex::Real /*x_width*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/) const;
};

struct TwoDimGaussianPulseFV
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*x*/,
        amrex::Real /*y*/,
        amrex::Real /*dx*/,
        amrex::Real /*dy*/,
        amrex::Real /*x0*/,
        amrex::Real /*y0*/,
        amrex::Real /*amplitude*/,
        amrex::Real /*x_width*/,
        amrex::Real /*y_width*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/) const;
};

struct SquarePulseFV
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*x*/,
        amrex::Real /*unused*/,
        amrex::Real /*dx*/,
        amrex::Real /*unused*/,
        amrex::Real /*x0*/,
        amrex::Real /*unused*/,
        amrex::Real /*amplitude*/,
        amrex::Real /*x_width*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/) const;
};

struct GaussianWavePacketFV
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*x*/,
        amrex::Real /*unused*/,
        amrex::Real /*dx*/,
        amrex::Real /*unused*/,
        amrex::Real /*x0*/,
        amrex::Real /*unused*/,
        amrex::Real /*amplitude*/,
        amrex::Real /*x_width*/,
        amrex::Real /*unused*/,
        amrex::Real /*x_wavenumber*/,
        amrex::Real /*unused*/) const;
};

} // namespace amr_wind

#endif /* ScalarAdvection_H */
