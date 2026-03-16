#ifndef W2AWAVES_H
#define W2AWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"
#include <vector>
#include <complex>
#ifdef AMR_WIND_USE_W2A
#include "Waves2AMR.h"
#include "AMReX_REAL.H"

using namespace amrex::literals;
#endif

namespace amr_wind::ocean_waves {

struct W2AWavesData : public RelaxZonesBaseData
{
    // Prefix for HOS files
    std::string modes_file{"modes_HOS_SWENSE.dat"};
    // Number of modes in each direction (related to FFT)
    int n0{0};
    int n1{0};
    int n2{0};
    // Number of points in each direction (spatial output)
    int n0_sp{0};
    int n1_sp{0};
    // Spatial resolution in each direction
    amrex::Real dx0{0.};
    amrex::Real dx1{0.};
    // Offsets for spatial data (between AMR-Wind domain and HOS domain)
    amrex::Real xlo0{0.};
    amrex::Real xlo1{0.};
    // Length to dimensionalize nondim data
    amrex::Real dimL{0.};
    // Timestep size of reference data
    amrex::Real dt_modes{0.0_rt};

    // Time (in HOS sim) to begin reading from
    // (this corresponds to t = 0 in AMR-Wind sim, functions as offset)
    amrex::Real t_winit{-1.0_rt};
    // Timestep (in HOS sim) of latest imported wave data
    int ntime{0};
    // Timestep (in HOS sim) of first imported wave data
    int n_winit{0};
    // Timestep (in HOS sim) of last available wave data
    int n_wstop{0};
    // Timestep offset (in HOS sim) for resetting to earlier data
    int n_offset{0};

    // Time (in AMR-Wind sim) when ow arrays were last changed
    // (this is from the last timestep, goes with time-interpolated data)
    amrex::Real t_last{-1.0_rt};
    // Time (in AMR-Wind sim) of latest imported wave data
    amrex::Real t{0.0_rt};

    // Vectors to store complex or real modes
    std::vector<std::complex<amrex::Real>> c_mX, c_mY, c_mZ, c_mFS;
    std::vector<amrex::Real> r_mX, r_mY, r_mZ, r_mFS, r_mAdd;

// Struct variables that have special types
#ifdef AMR_WIND_USE_W2A
    // FFTW plan vector
    std::vector<fftw_plan> plan_vector;
    // FFTW pointers to store modes ready for ifft
    fftw_complex *c_eta_mptr = nullptr, *c_u_mptr = nullptr,
                 *c_v_mptr = nullptr, *c_w_mptr = nullptr;
    // FFTW pointers to store modes ready for ifft
    amrex::Real *r_eta_mptr = nullptr, *r_u_mptr = nullptr, *r_v_mptr = nullptr,
                *r_w_mptr = nullptr, *au_mptr = nullptr, *av_mptr = nullptr,
                *aw_mptr = nullptr, *plan_out = nullptr;
    // ReadModes objects
    ReadModes<std::complex<amrex::Real>> c_rmodes;
    ReadModes<amrex::Real> r_rmodes;
#endif

    // Height vector (where velocity is sampled) stuff
    amrex::Vector<amrex::Real> hvec;
    // Index vector (where hvec overlaps with local boxes)
    amrex::Vector<int> indvec;
    // Flag indicating whether the source simulation used HOS-Ocean or HOS-NWT
    bool is_ocean{true};
    // Flag indicating interpolation should take place on this processor
    bool do_interp{true};
    // Flag indicating regrid has occurred since the last resizing of vecs
    bool resize_flag{false};
    // Vectors for spatial data from transformed modes
    amrex::Gpu::DeviceVector<amrex::Real> sp_eta_vec, sp_u_vec, sp_v_vec,
        sp_w_vec;

    ~W2AWavesData()
    {
#ifdef AMR_WIND_USE_W2A
        if (c_eta_mptr) {
            delete[] c_eta_mptr;
        }
        if (c_u_mptr) {
            delete[] c_u_mptr;
        }
        if (c_v_mptr) {
            delete[] c_v_mptr;
        }
        if (c_w_mptr) {
            delete[] c_w_mptr;
        }
        if (r_eta_mptr) {
            delete[] r_eta_mptr;
        }
        if (r_u_mptr) {
            delete[] r_u_mptr;
        }
        if (r_v_mptr) {
            delete[] r_v_mptr;
        }
        if (r_w_mptr) {
            delete[] r_w_mptr;
        }
        if (au_mptr) {
            delete[] au_mptr;
        }
        if (av_mptr) {
            delete[] av_mptr;
        }
        if (aw_mptr) {
            delete[] aw_mptr;
        }
        if (plan_out) {
            delete[] plan_out;
        }
        for (int n = 0; n < (int)plan_vector.size(); ++n) {
            fftw_destroy_plan(plan_vector[n]);
        }
        plan_vector.clear();
#endif
    }
};

struct W2AWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = W2AWavesData;
    using DataType = OceanWavesDataHolder<W2AWaves>;

    static std::string identifier() { return "W2AWaves"; }
};

} // namespace amr_wind::ocean_waves

#endif /* W2AWAVES_H */
