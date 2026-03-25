#ifndef ABLFIELDINIT_H
#define ABLFIELDINIT_H

#include "amr-wind/core/Field.H"

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

/** Initialize ABL fields
 */
class ABLFieldInit
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    ABLFieldInit();

    void initialize_from_inputfile(); // default
    void initialize_from_netcdf();    // general wind/temperature profiles

    void operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        const amrex::Array4<amrex::Real>& density,
        const amrex::Array4<amrex::Real>& temperature) const;

    /** Add temperature perturbations
     *
     *  This uses amrex::Random and, therefore, executes the loop on CPU and
     *  pushes the field to device.
     */
    void perturb_temperature(
        const int lev, const amrex::Geometry& geom, Field& temperature) const;

    //! Flag indicating whether temperature field needs perturbations
    bool add_temperature_perturbations() const { return m_perturb_theta; }

    //! Initialize TKE field
    void init_tke(const amrex::Geometry& geom, amrex::MultiFab& tke) const;

private:
    //! Initial velocity components
    amrex::Vector<amrex::Real> m_vel;

    //! File name for velocity forcing time table
    std::string m_vel_timetable;

    //! Initial wind speed
    amrex::Real m_vel_speed;

    //! Initial wind direction (rad)
    amrex::Real m_vel_dir;

    ///@{
    //! Lookup table containing temperature variation w.r.t. height above
    //! terrain
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;

    //! Adding option for wind heights similar to temperature heights
    //! Speed-up RANS calculation using 1-D profile for flat surface
    bool m_initial_wind_profile{false};
    bool m_terrain_aligned_profile{false};
    std::string m_terrain_file{"terrain.amrwind"};
    amrex::Gpu::DeviceVector<amrex::Real> m_xterrain;
    amrex::Gpu::DeviceVector<amrex::Real> m_yterrain;
    amrex::Gpu::DeviceVector<amrex::Real> m_zterrain;
    //! File name for 1-D data file
    std::string m_1d_rans;
    amrex::Vector<amrex::Real> m_wind_heights;
    amrex::Vector<amrex::Real> m_u_values;
    amrex::Vector<amrex::Real> m_v_values;
    amrex::Vector<amrex::Real> m_tke_values;
    ///@}

    // Device copies of the above arrays
    amrex::Gpu::DeviceVector<amrex::Real> m_thht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_thvv_d;

    amrex::Gpu::DeviceVector<amrex::Real> m_windht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_u_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_v_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_tke_d;
    //! Initial density field
    amrex::Real m_rho;

    //! Amplitude of fluctuations in x-direction
    amrex::Real m_deltaU{1.0};

    //! Amplitude of fluctuations in y-direction
    amrex::Real m_deltaV{1.0};

    //! Number of sinusoidal waves in x-direction
    amrex::Real m_Uperiods{4.0};

    //! Number of sinusoidal waves in y-direction
    amrex::Real m_Vperiods{4.0};

    //! Reference height for velocity perturbations
    amrex::Real m_ref_height{50.0};

    //! Reference IB height
  amrex::Real m_ib_height{0.0};
  
    //! Amplitude of temperature perturbations
    amrex::Real m_delta_t{0.8};

    //! Mean for Gaussian number generator
    amrex::Real m_theta_gauss_mean{0.0};

    //! Variance for Gaussian number generator
    amrex::Real m_theta_gauss_var{1.0};

    //! Cutoff height for temperature fluctuations
    amrex::Real m_theta_cutoff_height{1.0e16};

    //! Initial value for tke field
    amrex::Real m_tke_init{0.1};

    //! Multiplicative factor for init tke profile
    amrex::Real m_tke_init_factor{0.4};

    //! Cutoff height for init tke profile
    amrex::Real m_tke_cutoff_height{250.};

    //! Top velocity
    amrex::RealArray m_top_vel{20.0};

    //! Bottom velocity
    amrex::RealArray m_bottom_vel{4.0};

    //! Initial linear velocity profile boolean
    bool m_linear_profile{false};

    //! Perturb initial velocity field with sinusoidal fluctuations
    bool m_perturb_vel{true};

    //! Perturb temperature field with random fluctuations
    bool m_perturb_theta{false};

    //! Initialize tke profile non-constant
    bool m_tke_init_profile{false};

    //! Profiles of u, v, and theta are read in from a NetCDF File
    bool m_init_uvtheta_profile{false};
};

} // namespace amr_wind

#endif /* ABLFIELDINIT_H */
