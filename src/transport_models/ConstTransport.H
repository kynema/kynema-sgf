#ifndef CONSTTRANSPORT_H
#define CONSTTRANSPORT_H

#include "src/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::transport {

/** Fluid transport model with constant properties
 *  \ingroup transport
 */
class ConstTransport : public TransportModel::Register<ConstTransport>
{
public:
    static constexpr bool constant_properties = true;

    static std::string identifier() { return "ConstTransport"; }

    explicit ConstTransport(const CFDSim& sim) : m_repo(sim.repo())
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity", m_mu);
        pp.query("laminar_prandtl", m_Pr);
        pp.query("turbulent_prandtl", m_Prt);

        // Backwards compatibility
        amrex::ParmParse pp_boussinesq_buoyancy("BoussinesqBuoyancy");
        amrex::ParmParse pp_abl("ABL");
        if (pp.contains("thermal_expansion_coefficient")) {
            pp.get("thermal_expansion_coefficient", m_constant_beta);
            if (pp_boussinesq_buoyancy.contains("thermal_expansion_coeff")) {
                amrex::Print()
                    << "WARNING: BoussinesqBuoyancy.thermal_expansion_coeff "
                       "option has been deprecated in favor of "
                       "transport.thermal_expansion_coefficient. Ignoring the "
                       "BoussinesqBuoyancy option in favor of the transport "
                       "option."
                    << '\n';
            }
        } else if (pp_boussinesq_buoyancy.contains("thermal_expansion_coeff")) {
            amrex::Print()
                << "WARNING: BoussinesqBuoyancy.thermal_expansion_coeff option "
                   "has been deprecated in favor of "
                   "transport.thermal_expansion_coefficient. Please replace "
                   "this option."
                << '\n';
            pp_boussinesq_buoyancy.get(
                "thermal_expansion_coeff", m_constant_beta);
        }

        if (pp.contains("reference_temperature")) {
            pp.get("reference_temperature", m_reference_temperature);
            if (pp_boussinesq_buoyancy.contains("reference_temperature")) {
                amrex::Print()
                    << "WARNING: BoussinesqBuoyancy.reference_temperature "
                       "option has been deprecated in favor of "
                       "transport.reference_temperature. Ignoring the "
                       "BoussinesqBuoyancy option in favor of the transport "
                       "option."
                    << '\n';
            } else if (pp_abl.contains("reference_temperature")) {
                amrex::Print()
                    << "WARNING: ABL.reference_temperature "
                       "option has been deprecated in favor of "
                       "transport.reference_temperature. Ignoring the "
                       "ABL option in favor of the transport "
                       "option."
                    << '\n';
            }
        } else if (pp_boussinesq_buoyancy.contains("reference_temperature")) {
            amrex::Print()
                << "WARNING: BoussinesqBuoyancy.reference_temperature option "
                   "has been deprecated in favor of "
                   "transport.reference_temperature. Please replace "
                   "this option."
                << '\n';
            pp_boussinesq_buoyancy.get(
                "reference_temperature", m_reference_temperature);
        } else if (pp_abl.contains("reference_temperature")) {
            amrex::Print() << "WARNING: ABL.reference_temperature option "
                              "has been deprecated in favor of "
                              "transport.reference_temperature. Please replace "
                              "this option."
                           << '\n';
            pp_abl.get("reference_temperature", m_reference_temperature);
        }
    }

    ~ConstTransport() override = default;

    [[nodiscard]] amrex::Real viscosity() const { return m_mu; }

    [[nodiscard]] amrex::Real thermal_diffusivity() const
    {
        return m_mu / m_Pr;
    }

    [[nodiscard]] amrex::Real laminar_prandtl() const { return m_Pr; }

    [[nodiscard]] amrex::Real turbulent_prandtl() const { return m_Prt; }

    static amrex::Real laminar_schmidt(const std::string& scalar_name)
    {
        amrex::ParmParse pp("transport");
        const std::string key = scalar_name + "_laminar_schmidt";
        amrex::Real lam_schmidt = 1.0_rt;
        pp.query(key, lam_schmidt);
        return lam_schmidt;
    }

    static amrex::Real turbulent_schmidt(const std::string& scalar_name)
    {
        amrex::ParmParse pp("transport");
        const std::string key = scalar_name + "_turbulent_schmidt";
        amrex::Real turb_schmidt = 1.0_rt;
        pp.query(key, turb_schmidt);
        return turb_schmidt;
    }

    //! Return the dynamic visocity field
    std::unique_ptr<ScratchField> mu() override
    {
        auto mu = m_repo.create_scratch_field(1, m_ngrow);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*mu)(lev).setVal(m_mu);
        }
        return mu;
    }

    //! Return the thermal diffusivity field
    std::unique_ptr<ScratchField> alpha() override
    {
        auto alpha = mu();
        amrex::Real inv_Pr = 1.0_rt / m_Pr;
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*alpha)(lev).mult(inv_Pr);
        }
        return alpha;
    }

    std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& scalar_name) override
    {
        amrex::Real lam_schmidt = laminar_schmidt(scalar_name);

        amrex::Real inv_schmidt = 1.0_rt / lam_schmidt;
        auto diff = mu();
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*diff)(lev).mult(inv_schmidt);
        }

        return diff;
    }

    //! Return the thermal expansion coefficient
    [[nodiscard]] std::unique_ptr<ScratchField> beta() const override
    {
        auto beta = m_repo.create_scratch_field(1, m_ngrow);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            beta_fill(lev, (*beta)(lev));
        }
        return beta;
    }

    [[nodiscard]] amrex::Real reference_temperature() const override
    {
        return m_reference_temperature;
    }

    void beta_fill(const int lev, amrex::MultiFab& beta_mf) const override
    {
        if (m_constant_beta > 0.0_rt) {
            beta_mf.setVal(m_constant_beta);
        } else if (m_repo.field_exists("reference_temperature")) {
            const auto& temp0 = m_repo.get_field("reference_temperature");
            auto const& beta_arrs = beta_mf.arrays();
            auto const& temp0_arrs = temp0(lev).const_arrays();
            amrex::ParallelFor(
                beta_mf, amrex::IntVect(0), 1,
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                    beta_arrs[nbx](i, j, k) = 1.0_rt / temp0_arrs[nbx](i, j, k);
                });
            amrex::Gpu::streamSynchronize();
        } else {
            beta_mf.setVal(1.0_rt / m_reference_temperature);
        }

        if (m_repo.field_exists("vof")) {
            const auto& vof = m_repo.get_field("vof");
            auto const& beta_arrs = beta_mf.arrays();
            auto const& vof_arrs = vof(lev).const_arrays();
            amrex::ParallelFor(
                beta_mf, amrex::IntVect(0), 1,
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                    if (vof_arrs[nbx](i, j, k) > constants::TIGHT_TOL) {
                        beta_arrs[nbx](i, j, k) = 0.0_rt;
                    }
                });
            amrex::Gpu::streamSynchronize();
        }
    }

    void
    ref_theta_fill(const int lev, amrex::MultiFab& ref_theta_mf) const override
    {
        if (m_reference_temperature < 0.0_rt) {
            amrex::Abort("Reference temperature was not set");
        }

        if (m_repo.field_exists("reference_temperature")) {
            const auto& temp0 = m_repo.get_field("reference_temperature");
            amrex::MultiFab::Copy(ref_theta_mf, temp0(lev), 0, 0, 1, 0);
        } else {
            ref_theta_mf.setVal(m_reference_temperature);
        }
    }

    //! Return the reference temperature
    [[nodiscard]] std::unique_ptr<ScratchField> ref_theta() const override
    {
        if (m_reference_temperature < 0.0_rt) {
            amrex::Abort("Reference temperature was not set");
        }

        auto ref_theta = m_repo.create_scratch_field(1, m_ngrow);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            ref_theta_fill(lev, (*ref_theta)(lev));
        }
        return ref_theta;
    }

private:
    //! Reference to the field repository (for creating scratch fields)
    FieldRepo& m_repo;

    //! (Laminar) dynamic viscosity
    amrex::Real m_mu{1.0e-5_rt};

    //! Prandtl number
    amrex::Real m_Pr{1.0_rt};

    //! Turbulent Prandtl number
    amrex::Real m_Prt{1.0_rt};

    //! Constant thermal expansion coefficient
    amrex::Real m_constant_beta{0.0_rt};

    //! Reference temperature
    amrex::Real m_reference_temperature{-1.0_rt};
};

} // namespace kynema_sgf::transport

#endif /* CONSTTRANSPORT_H */
