#ifndef OVERSET_TYPES_H
#define OVERSET_TYPES_H

#include "AMReX_Gpu.H"

namespace amr_wind::tioga_iface {

/** Wrapper to handle host/device views of the data for TIOGA overset
 */
template <typename T>
struct AmrDualArray
{
    using ArrayType = amrex::Gpu::DeviceVector<T>;
    using HostArrayType = amrex::Vector<T>;

    //! Data on device
    ArrayType d_view;

    //! Data on host
    HostArrayType h_view;

    AmrDualArray() : d_view(), h_view() {}

    //! Construct arrays on host/device with a given length
    explicit AmrDualArray(const size_t len) : d_view(len), h_view(len) {}

    //! Resize arrays on host and device
    void resize(const size_t len)
    {
        d_view.resize(len);
        h_view.resize(len);
    }

    //! Return the size of the arrays
    [[nodiscard]] size_t size() const { return d_view.size(); }

    //! Copy data from host to device
    void copy_to_device()
    {
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_view.begin(), h_view.end(),
            d_view.begin());
    }

    //! Copy data from device to host
    void copy_to_host()
    {
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, d_view.begin(), d_view.end(),
            h_view.begin());
    }
};

} // namespace amr_wind::tioga_iface

#endif /* OVERSET_TYPES_H */
