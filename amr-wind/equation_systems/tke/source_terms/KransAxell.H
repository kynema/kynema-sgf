#ifndef KRANSAXELL_H
#define KRANSAXELL_H

#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/equation_systems/tke/TKESource.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::tke {

/** TKE source term based on Axell 2011 paper
 * Axell, L. B., & Liungman, O. (2001). A one-equation turbulence model for
 * geophysical applications: comparison with data and the k− ε model.
 * Environmental Fluid Mechanics, 1, 71-106.
 *  \ingroup tke_src turb_model we_abl
 */
class KransAxell : public TKESource::Register<KransAxell>
{
public:
    static std::string identifier() { return "KransAxell"; }

    explicit KransAxell(const CFDSim& /*sim*/);

    ~KransAxell() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_tke;
    amrex::Real m_Cmu{0.556_rt};
    amrex::Real m_heat_flux{0.0_rt};
    amrex::Real m_z0{0.1_rt};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_meso_start{600.0_rt};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    std::string m_1d_rans;
    amrex::Vector<amrex::Real> m_wind_heights;
    amrex::Vector<amrex::Real> m_tke_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_wind_heights_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_tke_values_d;
    //! Transport model
    const transport::TransportModel& m_transport;

    //! Reference temperature
    std::unique_ptr<ScratchField> m_ref_theta;
    bool m_horizontal_sponge{false};
    amrex::Real m_sponge_strength{1.0_rt};
    amrex::Real m_sponge_distance_west{-1000.0_rt};
    amrex::Real m_sponge_distance_east{1000.0_rt};
    amrex::Real m_sponge_distance_south{-1000.0_rt};
    amrex::Real m_sponge_distance_north{1000.0_rt};
    int m_sponge_west{0};
    int m_sponge_east{1};
    int m_sponge_south{0};
    int m_sponge_north{1};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
};

} // namespace amr_wind::pde::tke

#endif /* KRANSAXELL_H */
