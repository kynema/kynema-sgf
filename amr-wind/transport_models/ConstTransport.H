#ifndef CONSTTRANSPORT_H
#define CONSTTRANSPORT_H

#include "amr-wind/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::transport {

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
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi((*beta)(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& beta_arr = (*beta)(lev).array(mfi);
                beta_impl(lev, mfi, bx, beta_arr);
            }
        }
        return beta;
    }

    //! Compute the thermal expansion coefficient
    void beta_impl(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const amrex::Array4<amrex::Real>& beta) const override
    {

        if (m_constant_beta > 0.0_rt) {
            const amrex::Real beta_val = m_constant_beta;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                beta(i, j, k) = beta_val;
            });
        } else if (m_repo.field_exists("reference_temperature")) {
            const auto& temp0 = m_repo.get_field("reference_temperature");
            const auto& temp0_arr = temp0(lev).const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                beta(i, j, k) = 1.0_rt / temp0_arr(i, j, k);
            });
        } else {
            const amrex::Real beta_val = 1.0_rt / m_reference_temperature;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                beta(i, j, k) = beta_val;
            });
        }

        if (m_repo.field_exists("vof")) {
            const auto& vof = m_repo.get_field("vof");
            const auto& vof_arr = vof(lev).const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (vof_arr(i, j, k) > constants::TIGHT_TOL) {
                    beta(i, j, k) = 0.0_rt;
                }
            });
        }
    }

    [[nodiscard]] amrex::Real reference_temperature() const override
    {
        return m_reference_temperature;
    }

    //! Return the reference temperature
    [[nodiscard]] std::unique_ptr<ScratchField> ref_theta() const override
    {
        if (m_reference_temperature < 0.0_rt) {
            amrex::Abort("Reference temperature was not set");
        }

        auto ref_theta = m_repo.create_scratch_field(1, m_ngrow);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi((*ref_theta)(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& ref_theta_arr = (*ref_theta)(lev).array(mfi);
                ref_theta_impl(lev, mfi, bx, ref_theta_arr);
            }
        }
        return ref_theta;
    }

    //! Compute the reference temperature
    void ref_theta_impl(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const amrex::Array4<amrex::Real>& ref_theta) const override
    {
        if (m_reference_temperature < 0.0_rt) {
            amrex::Abort("Reference temperature was not set");
        }

        if (m_repo.field_exists("reference_temperature")) {
            auto& temp0 = m_repo.get_field("reference_temperature");
            const auto& temp0_arr = temp0(lev).const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                ref_theta(i, j, k) = temp0_arr(i, j, k);
            });
        } else {
            const amrex::Real ref_theta_val = m_reference_temperature;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                ref_theta(i, j, k) = ref_theta_val;
            });
        }
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

} // namespace amr_wind::transport

#endif /* CONSTTRANSPORT_H */
