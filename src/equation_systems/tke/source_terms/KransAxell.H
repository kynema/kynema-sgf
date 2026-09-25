#ifndef KRANSAXELL_H
#define KRANSAXELL_H

#include "src/transport_models/TransportModel.H"
#include "src/equation_systems/tke/TKESource.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

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
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

    /** Cap the shear production where the pressure-gradient sensor of
     *  KLAxellSeparation fires
     *
     *  Removes g max(0, P - C_P eps) from the source, with the ramp gate g of
     *  the sensor and C_P the production cap ratio of the model.
     *
     *  \param lev AMR level
     *  \param src_term Source term to modify
     */
    void production_cap(int lev, amrex::MultiFab& src_term) const;

    /** Cap the shear production with the time-relaxed gate of
     *  KLAxellSeparation
     *
     *  Same cap as production_cap, with the gate taken from the field
     *  separation_gate, used when KLAxellSeparation_coeffs.gate_relaxation_time
     *  is positive.
     *
     *  \param lev AMR level
     *  \param src_term Source term to modify
     */
    void production_cap_relaxed(int lev, amrex::MultiFab& src_term) const;

    /** Increase the dissipation where the pressure-gradient sensor of
     *  KLAxellSeparation fires
     *
     *  Removes g (c_d - 1) eps from the source, with the ramp gate g of the
     *  sensor and c_d the destruction boost factor of the model.
     *
     *  \param lev AMR level
     *  \param src_term Source term to modify
     */
    void destruction_boost(int lev, amrex::MultiFab& src_term) const;

    /** Increase the dissipation with the time-relaxed gate of
     *  KLAxellSeparation
     *
     *  Same boost as destruction_boost, with the gate taken from the field
     *  separation_gate.
     *
     *  \param lev AMR level
     *  \param src_term Source term to modify
     */
    void destruction_boost_relaxed(int lev, amrex::MultiFab& src_term) const;

    /** Treat the dissipation implicitly in the TKE diffusion solve
     *
     *  With the linearization eps = (Cmu^3 sqrt(k) / L) k, moves the implicit
     *  part of the dissipation (all of it for the implicit solve, half for
     *  Crank-Nicolson) from the explicit source to the diagonal of the
     *  diffusion solve through tke_lhs_src_term, as rho dt Cmu^3 sqrt(k) / L
     *  because the source is multiplied by the density before the solve. Not
     *  applied to the old-time source, which the Godunov scheme uses as
     *  forcing of the face states. Steady states are unchanged. The boundary,
     *  sponge and terrain terms, the production cap and the destruction boost
     *  stay explicit.
     *
     *  \param lev AMR level
     *  \param fstate Field state of the source evaluation
     *  \param src_term Source term to modify
     */
    void implicit_dissipation(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const;

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
    amrex::Real m_forcing_time_factor{5.0_rt};
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
    mutable std::unique_ptr<ScratchField> m_ref_theta_scratch;
    bool m_horizontal_sponge{false};
    amrex::Real m_sponge_strength{1.0_rt};
    amrex::Real m_sponge_distance_west{0.0_rt};
    amrex::Real m_sponge_distance_east{0.0_rt};
    amrex::Real m_sponge_distance_south{0.0_rt};
    amrex::Real m_sponge_distance_north{0.0_rt};
    bool m_sponge_west{false};
    bool m_sponge_east{false};
    bool m_sponge_south{false};
    bool m_sponge_north{false};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};

    //! Cap the shear production where the sensor of KLAxellSeparation fires
    bool m_production_cap{false};

    //! KLAxellSeparation stores a time-relaxed gate (separation_gate)
    bool m_relaxed_gate{false};

    //! Increase the dissipation where the sensor of KLAxellSeparation fires
    bool m_destruction_boost{false};

    //! Treat the dissipation implicitly (KLAxellSeparation)
    bool m_implicit_dissipation{false};

    //! Implicit fraction of the dissipation (1 implicit, 0.5 Crank-Nicolson)
    amrex::Real m_implicit_factor{1.0};
};

} // namespace kynema_sgf::pde::tke

#endif /* KRANSAXELL_H */
