#ifndef W2AWAVES_H
#define W2AWAVES_H

#include "src/ocean_waves/relaxation_zones/RelaxationZones.H"
#include <vector>
#include <complex>
#include <memory>
#ifdef KYNEMA_SGF_USE_W2A
#include "Waves2AMR.h"
#include "AMReX_REAL.H"

using namespace amrex::literals;
#endif

namespace kynema_sgf::ocean_waves {

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
    int max_padding_level{-1};
    // Spatial resolution in each direction
    amrex::Real dx0{0.};
    amrex::Real dx1{0.};
    // Offsets for spatial data (between Kynema-SGF domain and HOS domain)
    amrex::Real xlo0{0.};
    amrex::Real xlo1{0.};
    // Length to dimensionalize nondim data
    amrex::Real dimL{0.};
    // Timestep size of reference data
    amrex::Real dt_modes{0.0_rt};

    // Time (in HOS sim) to begin reading from
    // (this corresponds to t = 0 in Kynema-SGF sim, functions as offset)
    amrex::Real t_winit{-1.0_rt};
    // Timestep (in HOS sim) of latest imported wave data
    int ntime{0};
    // Timestep (in HOS sim) of first imported wave data
    int n_winit{0};
    // Timestep (in HOS sim) of last available wave data
    int n_wstop{0};
    // Timestep offset (in HOS sim) for resetting to earlier data
    int n_offset{0};

    // Time (in Kynema-SGF sim) of the last target update; negative before the
    // first update and zero after initialization with waves
    amrex::Real t_last{-1.0_rt};
    // Time (in Kynema-SGF sim) of latest imported wave data (New state of
    // w2a fields); the Old state holds the data at t - dt_modes
    amrex::Real t{0.0_rt};

    // Vectors to store complex or real modes
    std::vector<std::complex<amrex::Real>> c_mX, c_mY, c_mZ, c_mFS;
    std::vector<amrex::Real> r_mX, r_mY, r_mZ, r_mFS, r_mAdd;
    // Time-derivative modes for cubic Hermite time interpolation
    std::vector<std::complex<amrex::Real>> c_mT, c_mFST, c_mXt, c_mYt, c_mZt;
    std::vector<amrex::Real> r_mT, r_mFST, r_mAddT, r_mXt, r_mYt, r_mZt;
    // Cubic Hermite (with HOS time derivatives) or linear time interpolation
    bool cubic_time{true};

// Struct variables that have special types
#ifdef KYNEMA_SGF_USE_W2A
    // ReadModes objects
    ReadModes<std::complex<amrex::Real>> c_rmodes;
    ReadModes<amrex::Real> r_rmodes;
    modes_hosgrid::OceanReconstructionOptions reconstruction_options;
    interp_to_mfab::Interpolation interpolation{
        interp_to_mfab::Interpolation::cubic};
    std::unique_ptr<modes_hosgrid::OceanReconstruction> reconstruction;
    std::unique_ptr<modes_hosgrid::NwtReconstruction> nwt_reconstruction;
    std::size_t reported_selection{0};
    std::size_t reconstruction_reserved_bytes{0};

    int loaded_modes_step{std::numeric_limits<int>::min()};
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
};

struct W2AWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = W2AWavesData;
    using DataType = OceanWavesDataHolder<W2AWaves>;

    static std::string identifier() { return "W2AWaves"; }
};

} // namespace kynema_sgf::ocean_waves

#endif /* W2AWAVES_H */
