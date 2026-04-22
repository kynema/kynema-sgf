#ifndef SYNTHETICTURBULENCE_H
#define SYNTHETICTURBULENCE_H

#include <string>
#include <cmath>
#include <memory>

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

namespace synth_turb {

struct MeanVelOp
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real /* height */) const
    {
        return m_ref_vel;
    }

    const amrex::Real m_ref_vel;
};

class MeanProfile
{
public:
    explicit MeanProfile(amrex::Real ref_vel, int shear_dir = 2)
        : m_op{ref_vel}, m_ref_vel(ref_vel), m_shear_dir(shear_dir)
    {}

    virtual ~MeanProfile() = default;

    [[nodiscard]] amrex::Real reference_velocity() const { return m_ref_vel; }

    [[nodiscard]] int shear_dir() const { return m_shear_dir; }

    [[nodiscard]] MeanVelOp device_instance() const { return m_op; }

private:
    MeanVelOp m_op;

    const amrex::Real m_ref_vel;
    const int m_shear_dir;
};

} // namespace synth_turb

struct SynthTurbData
{
    // Dimensions of the box
    vs::VectorT<int> box_dims;

    // Length of the boxes in each direction
    vs::Vector box_len;

    vs::Vector dx;

    // Reference point for the turbulence box. Reference is the mid-point of the
    // turbulence grid at the plane where it is injected into the CFD flow
    // field.
    vs::Vector origin;

    // Transformation matrix to convert from global coordinate system to local
    // coordinate system
    vs::Tensor tr_mat;

    // Perturbation velocities (2, ny, nz)
    amrex::Vector<amrex::Real> uvel;
    amrex::Vector<amrex::Real> vvel;
    amrex::Vector<amrex::Real> wvel;

    amrex::Gpu::DeviceVector<amrex::Real> uvel_d;
    amrex::Gpu::DeviceVector<amrex::Real> vvel_d;
    amrex::Gpu::DeviceVector<amrex::Real> wvel_d;

    // Indices of the two planes stored in the data arrays
    int ileft;
    int iright;
};

struct SynthTurbDeviceData
{
    // Dimensions of the box
    vs::VectorT<int> box_dims;

    // Length of the boxes in each direction
    vs::Vector box_len;

    vs::Vector dx;

    // Reference point for the turbulence box. Reference is the mid-point of the
    // turbulence grid at the plane where it is injected into the CFD flow
    // field.
    vs::Vector origin;

    // Transformation matrix to convert from global coordinate system to local
    // coordinate system
    vs::Tensor tr_mat;

    amrex::Real* uvel;
    amrex::Real* vvel;
    amrex::Real* wvel;

    // Indices of the two planes stored in the data arrays
    int ileft;
    int iright;

    explicit SynthTurbDeviceData(SynthTurbData& hdata)
        : box_dims(hdata.box_dims)
        , box_len(hdata.box_len)
        , dx(hdata.dx)
        , origin(hdata.origin)
        , tr_mat(hdata.tr_mat)
        , uvel(hdata.uvel_d.data())
        , vvel(hdata.vvel_d.data())
        , wvel(hdata.wvel_d.data())
        , ileft(hdata.ileft)
        , iright(hdata.iright)
    {}
};

/** Indices and interpolation weights for a given point located within the
 *  turbulence box
 */
struct InterpWeights
{
    int il, ir, jl, jr, kl, kr;
    amrex::Real xl, xr, yl, yr, zl, zr;
};

/**
 *  \defgroup synth_turb Synthetic Turbulence
 *
 *  SyntheticTurbulence contains functions to inject turbulence into a
 *  CFD simulation. It reads a file to get the turbulent velocity on a
 *  grid and populates a source term for the momentum equation to
 *  achieve the desired turbulence characteristics at a specified
 *  plane in the CFD domain.
 *
 * \ingroup physics
 */

class SyntheticTurbulence : public Physics::Register<SyntheticTurbulence>
{
public:
    static std::string identifier() { return "SyntheticTurbulence"; }

    explicit SyntheticTurbulence(const CFDSim& sim);

    SyntheticTurbulence() = delete;
    SyntheticTurbulence(const SyntheticTurbulence&) = delete;
    SyntheticTurbulence& operator=(const SyntheticTurbulence&) = delete;

    ~SyntheticTurbulence() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override {}

    void initialize();

    void update();

    template <typename T>
    void update_impl(
        const SynthTurbDeviceData& /*turb_grid*/,
        const InterpWeights& /*weights*/,
        const T& /*velfunc*/);

private:
    const amr_wind::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_density;
    Field& m_turb_force;

    // Turbulence file name
    std::string m_turb_filename;

    // Turbulence box data
    SynthTurbData m_turb_grid;

    std::unique_ptr<synth_turb::MeanProfile> m_wind_profile;

    std::string m_mean_wind_type{"ConstValue"};
    amrex::Real m_epsilon;
    amrex::Real m_gauss_scaling;

    amrex::Real m_time_offset{0.0_rt};
    // Duration of the synthetic turbulence injection. Unused if < 0
    amrex::Real m_duration{-1};

    bool m_is_init{true};
};

} // namespace amr_wind

#endif /* SyntheticTurbulence_H */
