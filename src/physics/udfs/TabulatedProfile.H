#ifndef TABULATEDPROFILE_H
#define TABULATEDPROFILE_H

#include "src/utilities/linear_interpolation.H"

#include "AMReX_Geometry.H"
#include "AMReX_Gpu.H"
#include "AMReX_Orientation.H"
#include "AMReX_REAL.H"

#include <cstddef>

using namespace amrex::literals;

namespace kynema_sgf {

class Field;

namespace udf {

/** Inflow profile tabulated as a function of height in a text file
 *
 *  Reads a vertical profile from a whitespace-separated text file and imposes
 *  it as a Dirichlet condition on inflow boundaries. The file holds one row per
 *  height, with the first column the height and the remaining columns the
 *  field values, e.g. ``z u v T`` or ``z u v T tke``.
 *
 *  A different profile may be given on each face, so that the inflow faces
 *  follow the wind direction. Faces without a profile fall back to the constant
 *  value specified for that face.
 */
struct TabulatedProfile
{
    //! Number of boundary faces of the domain
    static constexpr int nfaces = AMREX_SPACEDIM * 2;

    struct DeviceOp
    {
        //! Heights of every face profile, concatenated
        const amrex::Real* z{nullptr};

        //! Field values of every face profile, concatenated and interleaved
        //! over the components of the field
        const amrex::Real* vals{nullptr};

        //! Index into z at which each face profile begins; negative when that
        //! face has no profile and the constant value is used instead
        amrex::GpuArray<int, nfaces> offset = {-1, -1, -1, -1, -1, -1};

        //! Number of heights in each face profile
        amrex::GpuArray<int, nfaces> npts = {0, 0, 0, 0, 0, 0};

        //! Constant value used on faces without a profile
        amrex::GpuArray<amrex::Real, nfaces * AMREX_SPACEDIM> constval = {
            0.0_rt};

        //! Ground height at each face, so that a profile tabulated above
        //! ground can be used on a boundary sitting on raised ground
        amrex::GpuArray<amrex::Real, nfaces> zoffset = {0.0_rt};

        //! Number of components of the field being filled
        int ncomp{1};

        //! Direction along which the profile is tabulated
        int idir{2};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& geom,
            const amrex::Real /*unused*/,
            amrex::Orientation ori,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            const int face = static_cast<int>(ori);
            const int fcomp = orig_comp + comp;

            amrex::Real val;
            if (offset[face] < 0) {
                val = constval[(face * AMREX_SPACEDIM) + fcomp];
            } else {
                const auto* problo = geom.ProbLo();
                const auto* dx = geom.CellSize();
                const auto zco =
                    problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);
                const auto zeff = zco - zoffset[face];

                const auto* zbegin = z + offset[face];
                const auto* zend = zbegin + npts[face];
                const auto* ybegin =
                    vals + (static_cast<std::ptrdiff_t>(ncomp) * offset[face]);
                val = interp::linear(zbegin, zend, ybegin, zeff, ncomp, fcomp);
            }
            field(iv[0], iv[1], iv[2], dcomp + comp) = val;
        }
    };

    using DeviceType = DeviceOp;

    static std::string identifier() { return "TabulatedProfile"; }

    explicit TabulatedProfile(const Field& fld);

    /** Return the device instance with its pointers bound to this object
     *
     *  The pointers cannot be set once and for all in the constructor: this
     *  object is copied into the fill patch operator, and each copy owns its
     *  own device buffers.
     */
    [[nodiscard]] DeviceType device_instance() const
    {
        DeviceType op = m_op;
        op.z = m_z_d.data();
        op.vals = m_vals_d.data();
        return op;
    }

    DeviceOp m_op;

    //! Device copy of the concatenated profile heights
    amrex::Gpu::DeviceVector<amrex::Real> m_z_d;

    //! Device copy of the concatenated profile values
    amrex::Gpu::DeviceVector<amrex::Real> m_vals_d;
};

} // namespace udf
} // namespace kynema_sgf

#endif /* TABULATEDPROFILE_H */
