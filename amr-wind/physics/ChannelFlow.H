#ifndef CHANNELFLOW_H
#define CHANNELFLOW_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/boundary_conditions/wall_models/WallFunction.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace {
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real analytical_smagorinsky_profile(
    const amrex::Real h,
    const amrex::Real Cs,
    const amrex::Real dx,
    const amrex::Real rho,
    const amrex::Real mu,
    const amrex::Real dpdx,
    const amrex::Real C0,
    const amrex::Real C1)
{
    const amrex::Real Cs2 = Cs * Cs;
    const amrex::Real dx2 = dx * dx;
    const amrex::Real Cs2dx2 = Cs2 * dx2;
    return (-mu / (rho * Cs2dx2) * h) +
           (Cs2dx2 / (3 * dpdx) *
            std::pow(
                ((2.0_rt / Cs2dx2 * dpdx * h) +
                 ((mu / (rho * Cs2dx2)) * (mu / (rho * Cs2dx2))) + C1),
                3.0_rt / 2.0_rt)) +
           C0;
}

} // namespace

namespace amr_wind::channel_flow {

/** Channel Flow physics
 *  \ingroup physics
 */
class ChannelFlow : public Physics::Register<ChannelFlow>
{
public:
    static std::string identifier() { return "ChannelFlow"; }

    explicit ChannelFlow(CFDSim& sim);

    ~ChannelFlow() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    //! Use Reichardt function to initialize velocity profile
    template <typename IndexSelector>
    void initialize_fields(
        int level,
        const amrex::Geometry& geom,
        const IndexSelector& idxOp,
        int n_idx);

    template <typename IndexSelector>
    amrex::Real compute_error(const IndexSelector& idxOp);

    amrex::Real compute_analytical_smagorinsky_error();

    void output_error();

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

private:
    //! CFD simulation controller instance
    CFDSim& m_sim;

    const amr_wind::SimTime& m_time;
    FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    WallFunction m_wall_func;

    //! Wall normal direction - Default y direction
    int m_norm_dir{1};

    //! initial density value
    amrex::Real m_rho{1.0_rt};

    //! viscosity
    amrex::Real m_mu{1.0_rt};

    //! Re_tau
    amrex::Real m_re_tau{1000.0_rt};

    //! u_tau
    amrex::Real m_utau{0.01_rt};

    //! y_tau
    amrex::Real m_ytau{0.001_rt};

    //! initial tke value
    amrex::Real m_tke0{0.0_rt};

    //! initial sdr value
    amrex::Real m_sdr0{1000.0_rt};

    //! flag for perturbations of the initial condition
    bool m_perturb_vel{false};

    //! perturbation period (y-direction)
    amrex::Real m_perturb_y_period{1.0_rt};

    //! perturbation period (z-direction)
    amrex::Real m_perturb_z_period{1.0_rt};

    //! perturbation factor (fraction of utau)
    amrex::Real m_perturb_fac{0.1_rt};

    //! Von-Karman constant
    amrex::Real m_kappa{0.41_rt};

    //! Turbulence model
    std::string m_turbulence_model;

    bool m_laminar{false};

    bool m_mesh_mapping{false};

    amrex::Real m_mean_vel;

    bool m_half{false};

    //! flag for analytical smagorinsky test
    bool m_analytical_smagorinsky_test{false};

    //! Analytical Smagorinsky first coefficient
    amrex::Real m_C0{0.0_rt};

    //! Analytical Smagorinsky second coefficient
    amrex::Real m_C1{0.0_rt};

    //! Body forcing (x direction)
    amrex::Real m_dpdx{0.0_rt};

    //! direction of mean velocity
    int m_mean_vel_dir{0};

    //! output precision
    const int m_w{18};

    //! error log file
    std::string m_output_fname{"channel_flow.log"};
};
} // namespace amr_wind::channel_flow

#endif /* ChannelFlow_H */
