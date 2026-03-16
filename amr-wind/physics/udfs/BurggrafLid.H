#ifndef BURGGRAF_LID_H
#define BURGGRAF_LID_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class Field;

namespace udf {

struct BurggrafLid
{
    struct DeviceOp
    {
        int idir{0};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real /*unused*/,
            amrex::Orientation /*unused*/,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            const auto* problo = geom.ProbLo();
            const auto* dx = geom.CellSize();
            const auto x = problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);
            amrex::Real val =
                16.0_rt * (std::pow(x, 4.0_rt) - 2.0_rt * std::pow(x, 3.0_rt) +
                           std::pow(x, 2.0_rt));
            amrex::RealArray uvec = {val};

            field(iv[0], iv[1], iv[2], dcomp + comp) = uvec[orig_comp + comp];
        }
    };
    using DeviceType = DeviceOp;

    static std::string identifier() { return "BurggrafLid"; }

    explicit BurggrafLid(const Field& fld);

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* BURGGRAF_LID_H */
