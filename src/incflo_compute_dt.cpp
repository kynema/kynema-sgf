#include "amr-wind/incflo.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

#include <cmath>
#include <limits>
#include "AMReX_REAL.H"

using namespace amrex::literals;

/** Estimate the new timestep for adaptive timestepping algorithm
 *
 * Compute new \f$\Delta t\f$ by using the formula derived in "A Boundary
 * Condition Capturing Method for Multiphase Incompressible Flow" by Kang et al.
 * (JCP).
 *
 *  \f{align}
 *  \text{CFL} &= \frac{\Delta t}{2}  \left[\left(C+V \right) +
 * \sqrt{\left(C+V\right)^2 + 4 \left(\frac{|F_x|}{\Delta x} +
 * \frac{|F_y|}{\Delta y} + \frac{|F_z|}{\Delta
 * z} \right)} \right] && \\
 *  C &= \text{max} \left(\frac{|U_x|}{\Delta x} ,\ \frac{|U_y|}{\Delta y} ,\
 * \frac{|U_z|}{\Delta z}\right) && \text{Convection} \\ V &= 2
 * \left(\frac{\mu}{\rho}\right)_\mathrm{max} \left[\frac{1}{\Delta x^2} +
 * \frac{1}{\Delta y^2} + \frac{1}{\Delta z^2} \right] && \text{Diffusion} \f}
 *
 *  \f$F_x, F_y, F_z\f$ are acceleration due to forcing terms. By default,
 *  \f$C\f$ is always used for computing CFL. The term \f$V\f$ is only used when
 *  the user has chosen implicit or Crank-Nicolson scheme for diffusion, and
 *  contributions from forcing term when `time.use_force_cfl` is `true` (default
 *  is `true`).
 *
 */
void incflo::compute_dt()
{
    BL_PROFILE("amr-wind::incflo::compute_dt");

    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);

    amrex::Real conv_cfl = 0.0_rt;
    amrex::Real diff_cfl = 0.0_rt;
    amrex::Real force_cfl = 0.0_rt;
    const bool mesh_mapping = m_sim.has_mesh_mapping();

    const auto& den = density();
    amr_wind::Field const* mesh_fac =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;
    const auto& mask_cell = m_repo.get_int_field("mask_cell");

    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const dxinv = geom[lev].InvCellSizeArray();
        amrex::MultiFab const& vel = icns().fields().field(lev);
        amrex::MultiFab const& vel_force = icns().fields().src_term(lev);
        amrex::MultiFab const& mu = icns().fields().mueff(lev);
        amrex::MultiFab const& rho = den(lev);

        auto const& vel_arr = vel.const_arrays();
        auto const& mask_arr = mask_cell(lev).const_arrays();
        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        amrex::Real conv_lev = 0.0_rt;
        amrex::Real mphase_conv_lev = 0.0_rt;
        amrex::Real diff_lev = 0.0_rt;
        amrex::Real force_lev = 0.0_rt;

        conv_lev += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real>{}, vel, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& v_bx = vel_arr[box_no];

                const amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0_rt;
                const amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0_rt;
                const amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0_rt;

                const auto mask =
                    static_cast<amrex::Real>(mask_arr[box_no](i, j, k));

                return amrex::max<amrex::Real>(
                    mask * std::abs(v_bx(i, j, k, 0)) * dxinv[0] / fac_x,
                    mask * std::abs(v_bx(i, j, k, 1)) * dxinv[1] / fac_y,
                    mask * std::abs(v_bx(i, j, k, 2)) * dxinv[2] / fac_z,
                    static_cast<amrex::Real>(-1.0_rt));
            });

        if (m_sim.pde_manager().has_pde("VOF")) {
            amrex::MultiFab const& vof = m_repo.get_field("vof")(lev);
            auto const& vof_arr = vof.const_arrays();
            mphase_conv_lev += amrex::ParReduce(
                amrex::TypeList<amrex::ReduceOpMax>{},
                amrex::TypeList<amrex::Real>{}, vel, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                    -> amrex::GpuTuple<amrex::Real> {
                    auto const& v_bx = vel_arr[box_no];
                    auto const& vof_bx = vof_arr[box_no];

                    // Check for interface
                    auto is_near =
                        amr_wind::multiphase::interface_band(i, j, k, vof_bx);

                    // CFL calculation is not needed away from interface
                    amrex::Real result = 0.0_rt;
                    if (is_near) {
                        // Near interface, evaluate CFL by sum of velocities
                        const amrex::Real fac_x =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 0))
                                         : 1.0_rt;
                        const amrex::Real fac_y =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 1))
                                         : 1.0_rt;
                        const amrex::Real fac_z =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 2))
                                         : 1.0_rt;

                        result =
                            (std::abs(v_bx(i, j, k, 0)) * dxinv[0] / fac_x) +
                            (std::abs(v_bx(i, j, k, 1)) * dxinv[1] / fac_y) +
                            (std::abs(v_bx(i, j, k, 2)) * dxinv[2] / fac_z);

                        const auto mask =
                            static_cast<amrex::Real>(mask_arr[box_no](i, j, k));

                        // Multiply advective CFL by 2 when near interface
                        result *= 2.0_rt * mask;
                        // CFL requirement to ensure vof conservation is 0.5_rt;
                        // this is half the typical concept of a CFL (1.0_rt)
                    }
                    return result;
                });
        }
        conv_lev = amrex::max(conv_lev, mphase_conv_lev);

        if (explicit_diffusion) {
            auto const& mu_arr = mu.const_arrays();
            auto const& rho_arr = rho.const_arrays();
            diff_lev += amrex::ParReduce(
                amrex::TypeList<amrex::ReduceOpMax>{},
                amrex::TypeList<amrex::Real>{}, rho, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                    -> amrex::GpuTuple<amrex::Real> {
                    auto const& mu_bx = mu_arr[box_no];
                    auto const& rho_bx = rho_arr[box_no];

                    const amrex::Real fac_x =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0_rt;
                    const amrex::Real fac_y =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0_rt;
                    const amrex::Real fac_z =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0_rt;

                    const amrex::Real dxinv2 =
                        2.0_rt * (dxinv[0] / fac_x * dxinv[0] / fac_x +
                                  dxinv[1] / fac_y * dxinv[1] / fac_y +
                                  dxinv[2] / fac_z * dxinv[2] / fac_z);

                    const auto mask =
                        static_cast<amrex::Real>(mask_arr[box_no](i, j, k));

                    return amrex::max<amrex::Real>(
                        mask * mu_bx(i, j, k) * dxinv2 / rho_bx(i, j, k),
                        -1.0_rt);
                });
        }

        if (m_time.use_force_cfl()) {
            auto const& vf_arr = vel_force.const_arrays();
            auto const& rho_arr = rho.const_arrays();
            force_lev += amrex::ParReduce(
                amrex::TypeList<amrex::ReduceOpMax>{},
                amrex::TypeList<amrex::Real>{}, vel_force, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                    -> amrex::GpuTuple<amrex::Real> {
                    auto const& vf_bx = vf_arr[box_no];
                    auto const& rho_bx = rho_arr[box_no];

                    const amrex::Real fac_x =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0_rt;
                    const amrex::Real fac_y =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0_rt;
                    const amrex::Real fac_z =
                        mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0_rt;

                    const auto mask =
                        static_cast<amrex::Real>(mask_arr[box_no](i, j, k));

                    return amrex::max<amrex::Real>(
                        mask * std::abs(vf_bx(i, j, k, 0)) * dxinv[0] / fac_x /
                            rho_bx(i, j, k),
                        mask * std::abs(vf_bx(i, j, k, 1)) * dxinv[1] / fac_y /
                            rho_bx(i, j, k),
                        mask * std::abs(vf_bx(i, j, k, 2)) * dxinv[2] / fac_z /
                            rho_bx(i, j, k),
                        static_cast<amrex::Real>(-1.0_rt));
                });
        }

        conv_cfl = amrex::max(conv_cfl, conv_lev);
        diff_cfl = amrex::max(diff_cfl, diff_lev);
        force_cfl = amrex::max(force_cfl, force_lev);
    }

    amrex::ParallelAllReduce::Max<amrex::Real>(
        conv_cfl, amrex::ParallelContext::CommunicatorSub());
    if (explicit_diffusion) {
        amrex::ParallelAllReduce::Max<amrex::Real>(
            diff_cfl, amrex::ParallelContext::CommunicatorSub());
    }
    if (m_time.use_force_cfl()) {
        amrex::ParallelAllReduce::Max<amrex::Real>(
            force_cfl, amrex::ParallelContext::CommunicatorSub());
    }

    m_time.set_current_cfl(conv_cfl, diff_cfl, force_cfl);
}

void incflo::compute_prescribe_dt()
{
    BL_PROFILE("amr-wind::incflo::compute_prescribe_dt");

    amrex::Real conv_cfl = 0.0_rt;
    const bool mesh_mapping = m_sim.has_mesh_mapping();

    amr_wind::Field const* mesh_fac =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;
    const auto& mask_cell = m_repo.get_int_field("mask_cell");

    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const dxinv = geom[lev].InvCellSizeArray();
        auto const& uf_arr = m_repo.get_field("u_mac")(lev).const_arrays();
        auto const& vf_arr = m_repo.get_field("v_mac")(lev).const_arrays();
        auto const& wf_arr = m_repo.get_field("w_mac")(lev).const_arrays();
        auto const& mask_arr = mask_cell(lev).const_arrays();

        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        amrex::Real conv_lev = 0.0_rt;
        amrex::Real mphase_conv_lev = 0.0_rt;

        conv_lev += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real>{}, icns().fields().field(lev),
            amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& umac = uf_arr[box_no];
                auto const& vmac = vf_arr[box_no];
                auto const& wmac = wf_arr[box_no];

                const amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0_rt;
                const amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0_rt;
                const amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0_rt;
                const auto mask =
                    static_cast<amrex::Real>(mask_arr[box_no](i, j, k));

                return amrex::max<amrex::Real>(
                    amrex::max<amrex::Real>(
                        std::abs(umac(i, j, k)), std::abs(umac(i + 1, j, k))) *
                        mask * dxinv[0] / fac_x,
                    amrex::max<amrex::Real>(
                        std::abs(vmac(i, j, k)), std::abs(vmac(i, j + 1, k))) *
                        mask * dxinv[1] / fac_y,
                    amrex::max<amrex::Real>(
                        std::abs(wmac(i, j, k)), std::abs(wmac(i, j, k + 1))) *
                        mask * dxinv[2] / fac_z,
                    static_cast<amrex::Real>(-1.0_rt));
            });

        if (m_sim.pde_manager().has_pde("VOF")) {
            amrex::MultiFab const& vof = m_repo.get_field("vof")(lev);
            auto const& vof_arr = vof.const_arrays();
            mphase_conv_lev += amrex::ParReduce(
                amrex::TypeList<amrex::ReduceOpMax>{},
                amrex::TypeList<amrex::Real>{}, icns().fields().field(lev),
                amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                    -> amrex::GpuTuple<amrex::Real> {
                    auto const& vof_bx = vof_arr[box_no];
                    auto const& umac = uf_arr[box_no];
                    auto const& vmac = vf_arr[box_no];
                    auto const& wmac = wf_arr[box_no];

                    // Check for interface
                    const auto is_near =
                        amr_wind::multiphase::interface_band(i, j, k, vof_bx);

                    // CFL calculation is not needed away from interface
                    amrex::Real result = 0.0_rt;
                    if (is_near) {
                        // Near interface, evaluate CFL by sum of velocities
                        const amrex::Real fac_x =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 0))
                                         : 1.0_rt;
                        const amrex::Real fac_y =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 1))
                                         : 1.0_rt;
                        const amrex::Real fac_z =
                            mesh_mapping ? (fac_arr[box_no](i, j, k, 2))
                                         : 1.0_rt;
                        const auto mask =
                            static_cast<amrex::Real>(mask_arr[box_no](i, j, k));

                        result = mask * (amrex::max(
                                             std::abs(umac(i, j, k)),
                                             std::abs(umac(i + 1, j, k))) *
                                             dxinv[0] / fac_x +
                                         amrex::max(
                                             std::abs(vmac(i, j, k)),
                                             std::abs(vmac(i, j + 1, k))) *
                                             dxinv[1] / fac_y +
                                         amrex::max(
                                             std::abs(wmac(i, j, k)),
                                             std::abs(wmac(i, j, k + 1))) *
                                             dxinv[2] / fac_z);
                    }
                    return result;
                });
        }
        conv_lev = amrex::max(conv_lev, mphase_conv_lev);

        conv_cfl = amrex::max(conv_cfl, conv_lev);
    }

    amrex::ParallelAllReduce::Max<amrex::Real>(
        conv_cfl, amrex::ParallelContext::CommunicatorSub());

    m_time.set_current_cfl(conv_cfl, 0.0_rt, 0.0_rt);
}
