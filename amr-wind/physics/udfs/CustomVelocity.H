#ifndef CUSTOM_VELOCITY_H
#define CUSTOM_VELOCITY_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class Field;

namespace udf {

struct CustomVelocity
{
    struct DeviceOp
    {
        // Declare parameters here if needed. For example:
        amrex::Real foo{1.0_rt};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> bar = {0.0_rt};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real /*time*/,
            amrex::Orientation /*ori*/,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            // Compute quantities to set the field values. For example:
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto x = problo[0] + ((iv[0] + 0.5_rt) * dx[0]);
            const auto y = problo[1] + ((iv[1] + 0.5_rt) * dx[1]);
            const auto z = problo[2] + ((iv[2] + 0.5_rt) * dx[2]);
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel = {
                1.0_rt, 0.0_rt, 0.0_rt};

            // Once the above is done, fill the field as:
            field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
            amrex::ignore_unused(x, y, z);
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "CustomVelocity"; }

    explicit CustomVelocity(const Field& fld);

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* CUSTOM_VELOCITY_H */
