#ifndef ABLISTATS_H
#define ABLISTATS_H

#include "amr-wind/wind_energy/ABLStatsBase.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "amr-wind/utilities/SecondMomentAveraging.H"
#include "amr-wind/utilities/ThirdMomentAveraging.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/wind_energy/ABLWallFunction.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

namespace pde::icns {
class ABLForcing;
class ABLMesoForcingMom;
} // namespace pde::icns

namespace pde::temperature {
class ABLMesoForcingTemp;
}

/**
 *  \defgroup abl_istats ABL Statistics
 *
 *  ABLStats contains functions to compute and statistics for ABL
 *  simulations. It supports output in ascii format as well as NetCDF format.
 *
 * \ingroup we_abl
 */
class ABLStats : public ABLStatsBase::Register<ABLStats>
{
public:
    static std::string identifier() { return "precursor"; }

    ABLStats(
        CFDSim& /*sim*/, const ABLWallFunction& /*abl_wall_func*/, int dir);

    ~ABLStats() override;

    ABLStatsMode abl_mode() const override { return ABLStatsMode::computed; }

    // Perform initialization actions after the mesh has been created
    void post_init_actions() override;

    //! Perform actions before a new timestep
    void pre_advance_work() override {}

    //! Process fields given timestep and output to disk
    void post_advance_work() override;

    //! Compute height of capping inversion
    void compute_zi();

    //! Return vel plane averaging instance
    const VelPlaneAveragingFine& vel_profile() const override
    {
        return m_pa_vel_fine;
    };

    //! Return vel plane averaging instance
    const VelPlaneAveraging& vel_profile_coarse() const override
    {
        return m_pa_vel;
    };

    //! Return instance that handles temperature statistics
    const FieldPlaneAveraging& theta_profile() const override
    {
        return m_pa_temp;
    }

    //! Return instance that handles temperature statistics
    const FieldPlaneAveragingFine& theta_profile_fine() const override
    {
        return m_pa_temp_fine;
    }

    void register_forcing_term(pde::icns::ABLForcing* forcing) const override
    {
        m_abl_forcing = forcing;
    }

    void register_meso_mom_forcing(
        pde::icns::ABLMesoForcingMom* forcing) const override
    {
        m_abl_meso_mom_forcing = forcing;
    }

    void register_meso_temp_forcing(
        pde::temperature::ABLMesoForcingTemp* forcing) const override
    {
        m_abl_meso_temp_forcing = forcing;
    }

    //! Calculate sfs stress averages
    void
    calc_sfs_stress_avgs(ScratchField& sfs_stress, ScratchField& t_sfs_stress);

    //! Calculate tke diffusion term by eliminating other terms
    void calc_tke_diffusion(
        ScratchField& diffusion,
        const Field& buoy_prod,
        const Field& shear_prod,
        const Field& dissipation,
        amrex::Real dt);

protected:
    //! Read user inputs and create the necessary files
    void initialize();

    //! Calculate plane average profiles
    void calc_averages();

    //! Output data based on user-defined format
    virtual void process_output();

    //! Prepare ASCII file
    virtual void prepare_ascii_file();
    //! Prepare NetCDF metadata
    virtual void prepare_netcdf_file();

    //! Write sampled data into a NetCDF file
    void write_netcdf();

    /** Output sampled data in ASCII format
     *
     *  Note that this should be used for debugging only and not in production
     *  runs as it can have significant impacts on code performance.
     */
    virtual void write_ascii();

private:
    CFDSim& m_sim;
    const ABLWallFunction& m_abl_wall_func;
    Field& m_temperature;
    Field& m_mueff;

    VelPlaneAveraging m_pa_vel;
    FieldPlaneAveraging m_pa_temp;
    VelPlaneAveragingFine m_pa_vel_fine;
    FieldPlaneAveragingFine m_pa_temp_fine;
    FieldPlaneAveraging m_pa_mueff;
    SecondMomentAveraging m_pa_tt;
    SecondMomentAveraging m_pa_tu;
    SecondMomentAveraging m_pa_uu;
    ThirdMomentAveraging m_pa_uuu;

    //! Reference to ABL forcing term if present
    mutable pde::icns::ABLForcing* m_abl_forcing{nullptr};

    mutable pde::icns::ABLMesoForcingMom* m_abl_meso_mom_forcing{nullptr};
    mutable pde::temperature::ABLMesoForcingTemp* m_abl_meso_temp_forcing{
        nullptr};

    //! Format of the data output (ascii, netcdf, etc.)
#ifdef AMR_WIND_USE_NETCDF
    std::string m_out_fmt{"netcdf"};
    std::string m_ncfile_name;
#else
    std::string m_out_fmt{"ascii"};
#endif
    std::string m_ascii_file_name;

    //! Frequency of data sampling and output
    int m_out_freq{100};

    //! Acceleration due to gravity magnitude
    amrex::Real m_gravity{9.81_rt};

    //! Von-Karman constant
    amrex::Real m_kappa{0.41_rt};

    //! Variable to store capping inversion height
    amrex::Real m_zi{0.0_rt};

    //! Wall-normal direction axis
    int m_normal_dir{2};

    //! Cell spacing at the coarsest level
    amrex::Real m_dn{0.0_rt};

    //! Number of cells in the horizontal direction
    size_t m_ncells_h1{0};
    size_t m_ncells_h2{0};

    //! Do energy budget
    bool m_do_energy_budget{false};
};

} // namespace amr_wind

#endif /* ABLISTATS_H */
