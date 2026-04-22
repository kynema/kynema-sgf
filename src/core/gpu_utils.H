#ifndef GPU_UTILS_H
#define GPU_UTILS_H

#include "AMReX_Gpu.H"

namespace kynema_sgf::gpu {

template <typename T>
amrex::Gpu::DeviceVector<T> device_view(const amrex::Vector<T>& hview)
{
    amrex::Gpu::DeviceVector<T> dview(hview.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, hview.begin(), hview.end(), dview.begin());
    return dview;
}

} // namespace kynema_sgf::gpu

#endif /* GPU_UTILS_H */
