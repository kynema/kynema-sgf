# Set amrex hydro options

set(HYDRO_SPACEDIM 3)
# Keep both flags set: HYDRO_EB=OFF skips EBGodunov/EBMOL targets, while
# HYDRO_NO_EB=ON disables EB code paths in AMReX-Hydro utility sources.
set(HYDRO_EB OFF CACHE BOOL "Disable AMReX-Hydro embedded boundary routines" FORCE)
set(HYDRO_NO_EB ON CACHE BOOL "Force-disable AMReX-Hydro EB code paths" FORCE)
set(HYDRO_MPI ${KYNEMA_SGF_ENABLE_MPI})
set(HYDRO_OMP ${KYNEMA_SGF_ENABLE_OPENMP})
set(HYDRO_FFT ${KYNEMA_SGF_ENABLE_FFT})
if (KYNEMA_SGF_ENABLE_CUDA)
  set(HYDRO_GPU_BACKEND CUDA CACHE STRING "HYDRO GPU type" FORCE)
elseif(KYNEMA_SGF_ENABLE_ROCM)
  set(HYDRO_GPU_BACKEND HIP CACHE STRING "HYDRO GPU type" FORCE)
elseif(KYNEMA_SGF_ENABLE_SYCL)
  set(HYDRO_GPU_BACKEND SYCL CACHE STRING "HYDRO GPU type" FORCE)
else()
  set(HYDRO_GPU_BACKEND NONE CACHE STRING "HYDRO GPU type" FORCE)
endif()
