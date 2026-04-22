#ifndef TWOPHASETRANSPORT_H
#define TWOPHASETRANSPORT_H

#include <numbers>
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::transport {

/** Two-phase fluid transport model with piece-wise constant properties
 *  \ingroup transport
 */
class TwoPhaseTransport : public TransportModel::Register<TwoPhaseTransport>
{
public:
    static constexpr bool constant_properties = false;

    static std::string identifier() { return "TwoPhaseTransport"; }

    explicit TwoPhaseTransport(const CFDSim& sim)
        : m_repo(sim.repo()), m_physics_mgr(sim.physics_manager())
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity_fluid1", m_mu1);
        pp.query("viscosity_fluid2", m_mu2);
        pp.query("laminar_prandtl_fluid1", m_Pr1);
        pp.query("laminar_prandtl_fluid2", m_Pr2);
        pp.query("turbulent_prandtl", m_Prt);

        // Check for single-phase quantities and warn
        if (pp.contains("viscosity")) {
            amrex::Print()
                << "WARNING: single-phase viscosity has been specified but "
                   "will not be used! (TwoPhaseTransport)\n";
        }
        if (pp.contains("laminar_prandtl")) {
            amrex::Print()
                << "WARNING: single-phase laminar_prandtl has been specified "
                   "but will not be used! (TwoPhaseTransport)\n";
        }

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

    ~TwoPhaseTransport() override = default;

    [[nodiscard]] amrex::Real laminar_prandtl1() const { return m_Pr1; }
    [[nodiscard]] amrex::Real laminar_prandtl2() const { return m_Pr2; }

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
        // Select the interface capturing method
        auto mu = m_repo.create_scratch_field(1, m_ngrow);

        auto ifacetype = get_iface_method();
        if (ifacetype == InterfaceCapturingMethod::VOF) {

            const auto& vof = m_repo.get_field("vof");

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                const auto& volfrac_arrs = vof(lev).const_arrays();
                const auto& visc_arrs = (*mu)(lev).arrays();
                const amrex::Real mu1 = m_mu1;
                const amrex::Real mu2 = m_mu2;
                amrex::ParallelFor(
                    (*mu)(lev), mu->num_grow(),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        visc_arrs[nbx](i, j, k) =
                            (mu1 * volfrac_arrs[nbx](i, j, k)) +
                            (mu2 * (1.0_rt - volfrac_arrs[nbx](i, j, k)));
                    });
            }
            amrex::Gpu::streamSynchronize();

        } else if (ifacetype == InterfaceCapturingMethod::LS) {

            const auto& levelset = m_repo.get_field("levelset");
            const auto& geom = m_repo.mesh().Geom();

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                const auto& dx = geom[lev].CellSizeArray();
                const auto& visc_arrs = (*mu)(lev).arrays();
                const auto& phi_arrs = levelset(lev).const_arrays();
                const amrex::Real eps =
                    std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);
                const amrex::Real mu1 = m_mu1;
                const amrex::Real mu2 = m_mu2;

                amrex::ParallelFor(
                    (*mu)(lev), mu->num_grow(),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        amrex::Real smooth_heaviside;
                        if (phi_arrs[nbx](i, j, k) > eps) {
                            smooth_heaviside = 1.0_rt;
                        } else if (phi_arrs[nbx](i, j, k) < -eps) {
                            smooth_heaviside = 0.0_rt;
                        } else {
                            smooth_heaviside =
                                0.5_rt *
                                (1.0_rt + phi_arrs[nbx](i, j, k) / eps +
                                 1.0_rt / std::numbers::pi_v<amrex::Real> *
                                     std::sin(
                                         phi_arrs[nbx](i, j, k) *
                                         std::numbers::pi_v<amrex::Real> /
                                         eps));
                        }
                        visc_arrs[nbx](i, j, k) =
                            (mu1 * smooth_heaviside) +
                            (mu2 * (1.0_rt - smooth_heaviside));
                    });
            }
            amrex::Gpu::streamSynchronize();
        }
        return mu;
    }

    //! Return the thermal diffusivity field (later divided by density, though)
    std::unique_ptr<ScratchField> alpha() override
    {
        // Select the interface capturing method
        auto alpha = m_repo.create_scratch_field(1, m_ngrow);

        auto ifacetype = get_iface_method();
        if (ifacetype == InterfaceCapturingMethod::VOF) {

            const auto& vof = m_repo.get_field("vof");

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                const auto& volfrac_arrs = vof(lev).const_arrays();
                const auto& thdiff_arrs = (*alpha)(lev).arrays();
                const amrex::Real mu1 = m_mu1;
                const amrex::Real mu2 = m_mu2;
                const amrex::Real Pr1 = m_Pr1;
                const amrex::Real Pr2 = m_Pr2;
                amrex::ParallelFor(
                    (*alpha)(lev), alpha->num_grow(),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        thdiff_arrs[nbx](i, j, k) =
                            (mu1 / Pr1 * volfrac_arrs[nbx](i, j, k)) +
                            (mu2 / Pr2 * (1.0_rt - volfrac_arrs[nbx](i, j, k)));
                    });
            }

        } else if (ifacetype == InterfaceCapturingMethod::LS) {
            const auto& levelset = m_repo.get_field("levelset");
            const auto& geom = m_repo.mesh().Geom();

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                const auto& dx = geom[lev].CellSizeArray();
                const auto& visc_arrs = (*alpha)(lev).arrays();
                const auto& phi_arrs = levelset(lev).const_arrays();
                const amrex::Real eps =
                    std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);
                const amrex::Real mu1 = m_mu1;
                const amrex::Real mu2 = m_mu2;
                const amrex::Real Pr1 = m_Pr1;
                const amrex::Real Pr2 = m_Pr2;
                amrex::ParallelFor(
                    (*alpha)(lev), alpha->num_grow(),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        amrex::Real smooth_heaviside;
                        if (phi_arrs[nbx](i, j, k) > eps) {
                            smooth_heaviside = 1.0_rt;
                        } else if (phi_arrs[nbx](i, j, k) < -eps) {
                            smooth_heaviside = 0.0_rt;
                        } else {
                            smooth_heaviside =
                                0.5_rt *
                                (1.0_rt + phi_arrs[nbx](i, j, k) / eps +
                                 1.0_rt / std::numbers::pi_v<amrex::Real> *
                                     std::sin(
                                         phi_arrs[nbx](i, j, k) *
                                         std::numbers::pi_v<amrex::Real> /
                                         eps));
                        }
                        visc_arrs[nbx](i, j, k) =
                            (mu1 / Pr1 * smooth_heaviside) +
                            (mu2 / Pr2 * (1.0_rt - smooth_heaviside));
                    });
            }
        }
        amrex::Gpu::streamSynchronize();
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

        auto ifacetype = get_iface_method();
        if (ifacetype == InterfaceCapturingMethod::VOF) {
            const auto& vof = m_repo.get_field("vof");
            const auto& vof_arr = vof(lev).const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (vof_arr(i, j, k) > constants::TIGHT_TOL) {
                    beta(i, j, k) = 0.0_rt;
                }
            });
        } else if (ifacetype == InterfaceCapturingMethod::LS) {
            const auto& levelset = m_repo.get_field("levelset");
            const auto& dx = m_repo.mesh().Geom(lev).CellSizeArray();
            const auto& phi_arr = levelset(lev).const_array(mfi);
            const amrex::Real eps = std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real smooth_heaviside;
                if (phi_arr(i, j, k) > eps) {
                    smooth_heaviside = 1.0_rt;
                } else if (phi_arr(i, j, k) < -eps) {
                    smooth_heaviside = 0.0_rt;
                } else {
                    smooth_heaviside =
                        0.5_rt *
                        (1.0_rt + phi_arr(i, j, k) / eps +
                         1.0_rt / std::numbers::pi_v<amrex::Real> *
                             std::sin(
                                 phi_arr(i, j, k) *
                                 std::numbers::pi_v<amrex::Real> / eps));
                }
                beta(i, j, k) *= (1.0_rt - smooth_heaviside);
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

    //! Reference to the physics manager
    const PhysicsMgr& m_physics_mgr;

    //! Phase 1 (liquid) dynamic molecular viscosity
    amrex::Real m_mu1{1.0e-3_rt};

    //! Phase 2 (gas) dynamic molecular viscosity
    amrex::Real m_mu2{1.0e-5_rt};

    //! Phase 1 (liquid) Prandtl number
    amrex::Real m_Pr1{7.2_rt};

    //! Phase 2 (gas) Prandtl number
    amrex::Real m_Pr2{0.7_rt};

    //! Turbulent Prandtl number
    amrex::Real m_Prt{1.0_rt};

    //! Constant thermal expansion coefficient
    amrex::Real m_constant_beta{0.0_rt};

    //! Reference temperature
    amrex::Real m_reference_temperature{-1.0_rt};

    [[nodiscard]] InterfaceCapturingMethod get_iface_method() const
    {
        if (!m_physics_mgr.contains("MultiPhase")) {
            amrex::Abort("TwoPhaseTransport requires MultiPhase physics");
        }
        const auto& multiphase = m_physics_mgr.get<MultiPhase>();
        return multiphase.interface_capturing_method();
    }
};

} // namespace amr_wind::transport

#endif /* TWOPHASETRANSPORT_H */
