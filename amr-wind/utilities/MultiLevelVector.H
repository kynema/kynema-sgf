#ifndef MULTILEVELVECTOR_H
#define MULTILEVELVECTOR_H

#include "AMReX_Geometry.H"
#include "AMReX_Vector.H"
#include "AMReX_GpuContainers.H"
#include "AMReX_MFParallelFor.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind {

/** Multi-level vector
 *  \ingroup utilities
 *
 *  This class enables multi-level vectors. The primary vector is the
 *  size of the number of levels and it contains host and device
 *  vectors sized for the number of points in a given direction of the
 *  domain.
 */
class MultiLevelVector
{
public:
    MultiLevelVector(const FieldLoc floc = FieldLoc::CELL) : m_floc(floc) {};

    void resize(int axis, const amrex::Vector<amrex::Geometry>& geom);

    [[nodiscard]] const amrex::Vector<amrex::Real>&
    host_data(const int lev) const
    {
        return m_data_h[lev];
    };

    amrex::Vector<amrex::Real>& host_data(const int lev)
    {
        return m_data_h[lev];
    };

    [[nodiscard]] const amrex::Gpu::DeviceVector<amrex::Real>&
    device_data(const int lev) const
    {
        return m_data_d[lev];
    };

    void copy_host_to_device();

    void copy_to_field(Field& fld);

    [[nodiscard]] int size() const { return static_cast<int>(m_data_h.size()); }

    [[nodiscard]] int ncells(const int lev) const
    {
        return static_cast<int>(m_data_h[lev].size());
    }

private:
    const FieldLoc m_floc;
    int m_axis{2};
    amrex::Vector<amrex::Gpu::DeviceVector<amrex::Real>> m_data_d;
    amrex::Vector<amrex::Vector<amrex::Real>> m_data_h;
};
} // namespace amr_wind
#endif
