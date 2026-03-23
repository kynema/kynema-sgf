#ifndef TWO_LAYER_H
#define TWO_LAYER_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class Field;

namespace udf {

struct TwoLayer
{
    struct DeviceOp
    {
        // velocities of the top and bottom layer respectively
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> top_vel = {0.0_rt};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> bottom_vel = {0.0_rt};
        amrex::Real z_part{0.5_rt}; // z-coordinate that divides top and bottom
        amrex::Real init_perturb{1.0_rt}; // initial perturbation for testing

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real time,
            amrex::Orientation /* ori */,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto z = problo[2] + ((iv[2] + 0.5_rt) * dx[2]);

            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel;
            if (z >= z_part) {
                vel = {top_vel[0], top_vel[1], 0.0_rt};
            } else {
                vel = {bottom_vel[0], bottom_vel[1], 0.0_rt};
            }

            if (time == 0.0_rt) {
                vel[0] *= init_perturb;
                vel[1] *= init_perturb;
            }

            field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "TwoLayer"; }

    explicit TwoLayer(const Field& fld);

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* TWO_LAYER_H */
