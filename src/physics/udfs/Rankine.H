#ifndef RANKINE_H
#define RANKINE_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

class Field;

namespace udf {

struct Rankine
{
    struct DeviceOp
    {
        amrex::Real Umax{15.0_rt};
        amrex::Real Rmax{125.0_rt};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> start_location = {
            AMREX_D_DECL(-10 * Rmax, 0.0_rt, 0.0_rt)};
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel_ref = {0.0_rt};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real time,
            amrex::Orientation /*unused*/,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto x = problo[0] + ((iv[0] + 0.5_rt) * dx[0]);
            const auto y = problo[1] + ((iv[1] + 0.5_rt) * dx[1]);
            const auto xr = x - (start_location[0] + vel_ref[0] * time);
            const auto yr = y - (start_location[1] + vel_ref[1] * time);
            const auto r = sqrt((xr * xr) + (yr * yr));
            const auto Vazi = Umax * ((r <= Rmax) ? (r / Rmax) : (Rmax / r));
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel = {
                AMREX_D_DECL(
                    static_cast<amrex::Real>(vel_ref[0] - Vazi * (yr / r)),
                    static_cast<amrex::Real>(vel_ref[1] + Vazi * (xr / r)),
                    static_cast<amrex::Real>(vel_ref[2] + 0.0_rt))};

            field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "Rankine"; }

    explicit Rankine(const Field& fld);

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace kynema_sgf

#endif /* RANKINE_H */
