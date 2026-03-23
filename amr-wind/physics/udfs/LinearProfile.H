#ifndef LINEARPROFILE_H
#define LINEARPROFILE_H

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class Field;

namespace udf {

struct LinearProfile
{
    struct DeviceOp
    {
        int idir{2};
        amrex::Real zmin;
        amrex::Real zmax;
        amrex::RealArray fld_min;
        amrex::RealArray fld_max;

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
            const auto zco = problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);

            amrex::Real val;
            if (zco < zmin) {
                val = fld_min[orig_comp + comp];
            } else if (zco > zmax) {
                val = fld_max[orig_comp + comp];
            } else {
                val = fld_min[orig_comp + comp] +
                      ((fld_max[orig_comp + comp] - fld_min[orig_comp + comp]) *
                       (zco - zmin) / (zmax - zmin));
            }
            field(iv[0], iv[1], iv[2], dcomp + comp) = val;
        }
    };

    using DeviceType = DeviceOp;

    static std::string identifier() { return "LinearProfile"; }

    explicit LinearProfile(const Field& fld);

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace udf
} // namespace amr_wind

#endif /* LINEARPROFILE_H */
