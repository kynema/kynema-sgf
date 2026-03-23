#include "amr-wind/physics/ChannelFlow.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/DirectionSelector.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::channel_flow {

ChannelFlow::ChannelFlow(CFDSim& sim)
    : m_sim(sim)
    , m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_wall_func(sim)
    , m_mesh_mapping(sim.has_mesh_mapping())
{
    {
        amrex::ParmParse pp("turbulence");
        pp.query("model", m_turbulence_model);
        if (m_turbulence_model == "Laminar") {
            m_laminar = true;
        }
    }
    {
        {
            amrex::ParmParse pp("incflo");
            pp.query("density", m_rho);
        }
        amrex::ParmParse pp("ChannelFlow");
        pp.query("flow_direction", m_mean_vel_dir);
        pp.query("normal_direction", m_norm_dir);
        pp.query("error_log_file", m_output_fname);

        if (m_laminar) {
            pp.query("Mean_Velocity", m_mean_vel);
            pp.query("half_channel", m_half);
        } else {
            pp.query("re_tau", m_re_tau);
            pp.query("tke0", m_tke0);
            pp.query("sdr0", m_sdr0);
            pp.query("perturb_velocity", m_perturb_vel);
            pp.query("perturb_y_period", m_perturb_y_period);
            pp.query("perturb_z_period", m_perturb_z_period);
            pp.query("perturb_factor", m_perturb_fac);
            pp.query(
                "analytical_smagorinsky_test", m_analytical_smagorinsky_test);
            if (m_analytical_smagorinsky_test) {
                if (m_turbulence_model == "Smagorinsky") {
                    pp.query("C0", m_C0);
                    pp.query("C1", m_C1);
                    {
                        amrex::ParmParse ppb("BodyForce");
                        amrex::Vector<amrex::Real> body_force{
                            0.0_rt, 0.0_rt, 0.0_rt};
                        ppb.queryarr("magnitude", body_force);
                        m_dpdx = -body_force[0];
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                            std::abs(body_force[1]) <
                                std::numeric_limits<amrex::Real>::epsilon(),
                            "body force in y should be zero for this wall "
                            "function");
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                            std::abs(body_force[2]) <
                                std::numeric_limits<amrex::Real>::epsilon(),
                            "body force in z should be zero for this wall "
                            "function");
                    }
                } else {
                    amrex::Abort(
                        "Use Smagorinsky model for the analytical Smagorinsky "
                        "test");
                }
            }
        }
    }
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity", m_mu);
        // Assumes a boundary layer height of 1.0_rt
        m_utau = m_mu * m_re_tau / (m_rho * 1.0_rt);
        m_ytau = m_mu / (m_utau * m_rho);
    }
    if ((amrex::ParallelDescriptor::IOProcessor()) &&
        (m_laminar || m_analytical_smagorinsky_test)) {
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u" << '\n';
        f.close();
    }

    if (!m_repo.field_exists("wall_dist")) {
        m_repo.declare_field("wall_dist", 1, 1, 1);
    }
}

/** Initialize the velocity, density, tke and sdr fields at the beginning of the
 *  simulation.
 */
void ChannelFlow::initialize_fields(int level, const amrex::Geometry& geom)
{
    switch (m_norm_dir) {
    case 0:
        initialize_fields(level, geom, XDir(), 0);
        break;
    case 1:
        initialize_fields(level, geom, YDir(), 1);
        break;
    case 2:
        initialize_fields(level, geom, ZDir(), 2);
        break;
    default:
        amrex::Abort("axis must be equal to 1 or 2");
        break;
    }
}

template <typename IndexSelector>
void ChannelFlow::initialize_fields(
    int level,
    const amrex::Geometry& geom,
    const IndexSelector& idxOp,
    const int n_idx)
{
    const amrex::Real kappa = m_kappa;
    const amrex::Real y_tau = m_ytau;
    const amrex::Real utau = m_utau;
    auto& velocity_field = m_repo.get_field("velocity");
    auto& velocity = velocity_field(level);
    auto& density = m_repo.get_field("density")(level);
    density.setVal(m_rho);
    const auto perturb_vel = m_perturb_vel;
    const auto perturb_fac = m_perturb_fac;
    const auto perturb_amp = perturb_fac * m_utau;
    const auto pi = std::numbers::pi_v<amrex::Real>;
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const auto& dx = geom.CellSizeArray();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> lengths{AMREX_D_DECL(
        probhi[0] - problo[0], probhi[1] - problo[1], probhi[2] - problo[2])};
    const auto y_perturb = m_perturb_y_period * 2.0_rt * pi / lengths[1];
    const auto z_perturb = m_perturb_z_period * 2.0_rt * pi / lengths[2];

    if (!m_laminar) {
        if (m_analytical_smagorinsky_test) {
            auto coeffs = m_sim.turbulence_model().model_coeffs();
            const auto Cs = coeffs["Cs"];
            const auto mu = m_mu;
            const auto rho = m_rho;
            const auto C0 = m_C0;
            const auto C1 = m_C1;
            const auto dpdx = m_dpdx;
            const auto& vel_arrs = velocity.arrays();

            amrex::ParallelFor(
                velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const int n_ind = idxOp(i, j, k);
                    amrex::Real h =
                        problo[n_idx] + ((n_ind + 0.5_rt) * dx[n_idx]);
                    if (h > 1.0_rt) {
                        h = 2.0_rt - h;
                    }
                    const amrex::Real ux = analytical_smagorinsky_profile(
                        h, Cs, dx[n_idx], rho, mu, dpdx, C0, C1);
                    vel_arrs[nbx](i, j, k, 0) =
                        ux + (perturb_vel
                                  ? perturb_fac * std::sin(z_perturb * h)
                                  : 0.0_rt);
                    vel_arrs[nbx](i, j, k, 1) = 0.0_rt;
                    vel_arrs[nbx](i, j, k, 2) = 0.0_rt;
                });
            amrex::Gpu::streamSynchronize();
        } else {
            if (m_repo.field_exists("tke")) {
                auto& tke = m_repo.get_field("tke")(level);
                tke.setVal(m_tke0);
            }
            if (m_repo.field_exists("sdr")) {
                auto& sdr = m_repo.get_field("sdr")(level);
                sdr.setVal(m_sdr0);
            }
            auto& walldist = m_repo.get_field("wall_dist")(level);

            const auto& vel_arrs = velocity.arrays();
            const auto& wd_arrs = walldist.arrays();

            amrex::ParallelFor(
                velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const int n_ind = idxOp(i, j, k);
                    amrex::Real h =
                        problo[n_idx] + ((n_ind + 0.5_rt) * dx[n_idx]);
                    if (h > 1.0_rt) {
                        h = 2.0_rt - h;
                    }
                    wd_arrs[nbx](i, j, k) = h;
                    const amrex::Real hp = h / y_tau;
                    vel_arrs[nbx](i, j, k, 0) =
                        utau *
                        (1.0_rt / kappa * std::log1p(kappa * hp) +
                         7.8_rt * (1.0_rt - std::exp(-hp / 11.0_rt) -
                                   (hp / 11.0_rt) * std::exp(-hp / 3.0_rt)));

                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                    const amrex::Real perty = z_perturb * perturb_amp *
                                              std::sin(y_perturb * y) *
                                              std::cos(z_perturb * z);
                    const amrex::Real pertz = -y_perturb * perturb_amp *
                                              std::cos(y_perturb * y) *
                                              std::sin(z_perturb * z);
                    vel_arrs[nbx](i, j, k, 1) =
                        0.0_rt + (perturb_vel ? perty : 0.0_rt);
                    vel_arrs[nbx](i, j, k, 2) =
                        0.0_rt + (perturb_vel ? pertz : 0.0_rt);
                });
            amrex::Gpu::streamSynchronize();
        }
    } else {
        velocity.setVal(0.0_rt);
        velocity.setVal(
            m_mean_vel, m_mean_vel_dir, 1,
            velocity_field.num_grow()[m_mean_vel_dir]);
    }
}

template <typename IndexSelector>
amrex::Real ChannelFlow::compute_error(const IndexSelector& idxOp)
{
    amrex::Real error = 0.0_rt;
    const auto flow_dir = m_mean_vel_dir;
    const auto norm_dir = m_norm_dir;
    const auto mu = m_mu;
    const auto mesh_mapping = m_mesh_mapping;

    amrex::Real dpdx = 0;
    {
        amrex::Vector<amrex::Real> body_force{0.0_rt, 0.0_rt, 0.0_rt};
        amrex::ParmParse pp("BodyForce");
        pp.queryarr("magnitude", body_force, 0, AMREX_SPACEDIM);
        dpdx = body_force[flow_dir];
    }

    const auto& problo = m_mesh.Geom(0).ProbLoArray();
    amrex::Real ht = 0.0_rt;
    {
        amrex::Vector<amrex::Real> probhi_physical{0.0_rt, 0.0_rt, 0.0_rt};
        amrex::ParmParse pp("geometry");
        if (pp.contains("prob_hi_physical")) {
            pp.getarr("prob_hi_physical", probhi_physical);
        } else {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                probhi_physical[d] = m_mesh.Geom(0).ProbHiArray()[d];
            }
        }
        ht = probhi_physical[norm_dir] - problo[norm_dir];
        // For half channel, channel height is double the domain height
        if (m_half) {
            ht *= 2.0_rt;
        }
    }

    const auto& velocity = m_repo.get_field("velocity");
    Field const* nu_coord_cc =
        mesh_mapping ? &(m_repo.get_field("non_uniform_coord_cc")) : nullptr;
    Field const* mesh_fac_cc =
        mesh_mapping
            ? &(m_repo.get_mesh_mapping_field(amr_wind::FieldLoc::CELL))
            : nullptr;

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& prob_lo = m_mesh.Geom(lev).ProbLoArray();

        const auto& vel = velocity(lev);
        auto const& vel_arr = vel.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();
        amrex::MultiArray4<amrex::Real const> fac_arr =
            mesh_mapping ? ((*mesh_fac_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();
        amrex::MultiArray4<amrex::Real const> nu_cc =
            mesh_mapping ? ((*nu_coord_cc)(lev).const_arrays())
                         : amrex::MultiArray4<amrex::Real const>();

        error += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real>{}, vel, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& vel_bx = vel_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                amrex::Real y =
                    mesh_mapping ? (nu_cc[box_no](i, j, k, norm_dir))
                                 : (prob_lo[norm_dir] +
                                    ((idxOp(i, j, k) + 0.5_rt) * dx[norm_dir]));
                amrex::Real fac_x =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 0)) : 1.0_rt;
                amrex::Real fac_y =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 1)) : 1.0_rt;
                amrex::Real fac_z =
                    mesh_mapping ? (fac_arr[box_no](i, j, k, 2)) : 1.0_rt;

                const amrex::Real u = vel_bx(i, j, k, flow_dir);
                const amrex::Real u_exact =
                    1.0_rt / (2.0_rt * mu) * -dpdx * (y * y - y * ht);

                const amrex::Real cell_vol =
                    dx[0] * fac_x * dx[1] * fac_y * dx[2] * fac_z;

                return cell_vol * mask_bx(i, j, k) * (u - u_exact) *
                       (u - u_exact);
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

amrex::Real ChannelFlow::compute_analytical_smagorinsky_error()
{
    amrex::Real error = 0.0_rt;
    const auto mu = m_mu;
    auto coeffs = m_sim.turbulence_model().model_coeffs();
    const auto Cs = coeffs["Cs"];
    const auto rho = m_rho;
    const auto C0 = m_C0;
    const auto C1 = m_C1;
    const auto dpdx = m_dpdx;

    const auto& velocity = m_repo.get_field("velocity");

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& prob_lo = m_mesh.Geom(lev).ProbLoArray();

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& vel = velocity(lev);
        auto const& vel_arr = vel.const_arrays();
        auto const& mask_arr = level_mask.const_arrays();

        error += amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real>{}, vel, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                auto const& vel_bx = vel_arr[box_no];
                auto const& mask_bx = mask_arr[box_no];

                const int n_idx = 2;
                amrex::Real h = prob_lo[n_idx] + ((k + 0.5_rt) * dx[n_idx]);
                if (h > 1.0_rt) {
                    h = 2.0_rt - h;
                }
                const amrex::Real u = vel_bx(i, j, k, 0);
                const amrex::Real u_exact = analytical_smagorinsky_profile(
                    h, Cs, dx[n_idx], rho, mu, dpdx, C0, C1);

                const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

                return cell_vol * mask_bx(i, j, k) * (u - u_exact) *
                       (u - u_exact);
            });
    }

    amrex::ParallelDescriptor::ReduceRealSum(error);

    const auto total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void ChannelFlow::output_error()
{
    // analytical solution exists only for flow in streamwise direction
    amrex::Real u_err = 0.0_rt;
    switch (m_norm_dir) {
    case 0:
        u_err = compute_error(XDir());
        break;
    case 1:
        u_err = compute_error(YDir());
        break;
    case 2:
        if (m_laminar) {
            u_err = compute_error(ZDir());
        } else if (m_analytical_smagorinsky_test) {
            u_err = compute_analytical_smagorinsky_error();
        }
        break;
    default:
        amrex::Abort("axis must be equal to 1 or 2");
        break;
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.new_time()
          << std::setw(m_w) << u_err << '\n';
        f.close();
    }
}

void ChannelFlow::post_init_actions()
{
    if (m_laminar || m_analytical_smagorinsky_test) {
        output_error();
    }

    auto& velocity = m_repo.get_field("velocity");
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if ((velocity.bc_type()[zlo] == BC::wall_model) ||
        (velocity.bc_type()[zhi] == BC::wall_model)) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_norm_dir == 2,
            "Wall normal direction should be 2 if using a wall model");
        velocity.register_custom_bc<VelWallFunc>(m_wall_func);

        amrex::ParmParse pp("BodyForce");
        amrex::Vector<amrex::Real> body_force{0.0_rt, 0.0_rt, 0.0_rt};
        pp.getarr("magnitude", body_force);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs(body_force[2]) <
                std::numeric_limits<amrex::Real>::epsilon(),
            "body force in z should be zero if using a wall model");
    }
}

void ChannelFlow::post_advance_work()
{
    if (m_laminar || m_analytical_smagorinsky_test) {
        output_error();
    }
}

} // namespace amr_wind::channel_flow
