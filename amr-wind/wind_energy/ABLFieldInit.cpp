#include <cmath>
#include <numbers>
#include <string>
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX_ParallelDescriptor.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/io_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

ABLFieldInit::ABLFieldInit()
{
    amrex::ParmParse pp_abl("ABL");
    if (pp_abl.contains("init_profile")) {
        initialize_from_netcdf();
    } else {
        initialize_from_inputfile();
    }
}

void ABLFieldInit::initialize_from_inputfile()
{
    amrex::ParmParse pp_abl("ABL");
    //! Check for wind profile
    pp_abl.query("initial_wind_profile", m_initial_wind_profile);
    pp_abl.query("terrain_aligned_profile", m_terrain_aligned_profile);
    amrex::ParmParse pp_terrain("TerrainDrag");
    pp_terrain.query("terrain_file", m_terrain_file);
    if (m_initial_wind_profile) {
        pp_abl.query("rans_1dprofile_file", m_1d_rans);
        if (!m_1d_rans.empty()) {
            std::ifstream ransfile(m_1d_rans, std::ios::in);
            if (!ransfile.good()) {
                amrex::Abort("Cannot find 1-D RANS profile file " + m_1d_rans);
            }
            amrex::Real value1, value2, value3, value4, value5;
            while (ransfile >> value1 >> value2 >> value3 >> value4 >> value5) {
                m_wind_heights.push_back(value1);
                m_u_values.push_back(value2);
                m_v_values.push_back(value3);
                m_tke_values.push_back(value5);
            }
        }
        amrex::Vector<amrex::Real> xterrain;
        amrex::Vector<amrex::Real> yterrain;
        amrex::Vector<amrex::Real> zterrain;
        if (m_terrain_aligned_profile) {
            ioutils::read_flat_grid_file(
                m_terrain_file, xterrain, yterrain, zterrain);
        }
        const auto xterrain_size = xterrain.size();
        const auto yterrain_size = yterrain.size();
        const auto zterrain_size = zterrain.size();
        m_xterrain.resize(xterrain_size);
        m_yterrain.resize(yterrain_size);
        m_zterrain.resize(zterrain_size);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, xterrain.begin(), xterrain.end(),
            m_xterrain.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, yterrain.begin(), yterrain.end(),
            m_yterrain.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, zterrain.begin(), zterrain.end(),
            m_zterrain.begin());
    }
    // Temperature variation as a function of height
    pp_abl.getarr("temperature_heights", m_theta_heights);
    pp_abl.getarr("temperature_values", m_theta_values);

    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    int num_theta_values = static_cast<int>(m_theta_heights.size());

    pp_abl.query("perturb_velocity", m_perturb_vel);
    pp_abl.query("perturb_ref_height", m_ref_height);
    pp_abl.query("Uperiods", m_Uperiods);
    pp_abl.query("Vperiods", m_Vperiods);
    pp_abl.query("deltaU", m_deltaU);
    pp_abl.query("deltaV", m_deltaV);

    pp_abl.query("perturb_temperature", m_perturb_theta);
    pp_abl.query("random_gauss_mean", m_theta_gauss_mean);
    pp_abl.query("random_gauss_var", m_theta_gauss_var);
    pp_abl.query("cutoff_height", m_theta_cutoff_height);
    pp_abl.query("theta_amplitude", m_delta_t);

    pp_abl.query("init_tke", m_tke_init);
    pp_abl.query("init_tke_beare_profile", m_tke_init_profile);
    pp_abl.query("init_tke_beare_factor", m_tke_init_factor);
    pp_abl.query("init_tke_cutoff_height", m_tke_cutoff_height);

    pp_abl.query("linear_profile", m_linear_profile);

    pp_abl.query("top_velocity", m_top_vel);
    pp_abl.query("bottom_velocity", m_bottom_vel);

    // TODO: Modify this to accept velocity as a function of height
    amrex::ParmParse pp_incflo("incflo");
    amrex::ParmParse pp_mphase("MultiPhase");
    if (!pp_mphase.contains("density_fluid2")) {
        pp_incflo.get("density", m_rho);
    } else {
        pp_mphase.get("density_fluid2", m_rho);
        // Note: density field will later be overwritten by MultiPhase
        // post_init
    }

    amrex::ParmParse pp_forcing("ABLForcing");
    pp_forcing.query("velocity_timetable", m_vel_timetable);
    if (!m_vel_timetable.empty()) {
        std::ifstream ifh(m_vel_timetable, std::ios::in);
        if (!ifh.good()) {
            amrex::Abort(
                "Cannot find ABLForcing velocity_timetable file: " +
                m_vel_timetable);
        }
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        amrex::Real vel_time;
        amrex::Real vel_ang;
        ifh >> vel_time >> m_vel_speed >> vel_ang;
        m_vel_dir = utils::radians(vel_ang);
    } else {
        pp_incflo.getarr("velocity", m_vel);
    }

    m_thht_d.resize(num_theta_values);
    m_thvv_d.resize(num_theta_values);
    if (m_initial_wind_profile) {
        int num_wind_values = static_cast<int>(m_wind_heights.size());
        m_windht_d.resize(num_wind_values);
        m_prof_u_d.resize(num_wind_values);
        m_prof_v_d.resize(num_wind_values);
        m_prof_tke_d.resize(num_wind_values);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_wind_heights.begin(),
            m_wind_heights.end(), m_windht_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_u_values.begin(), m_u_values.end(),
            m_prof_u_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_v_values.begin(), m_v_values.end(),
            m_prof_v_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_tke_values.begin(), m_tke_values.end(),
            m_prof_tke_d.begin());
    }
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_heights.begin(),
        m_theta_heights.end(), m_thht_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
        m_thvv_d.begin());
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void ABLFieldInit::initialize_from_netcdf()
{
#ifndef AMR_WIND_USE_NETCDF
    amrex::Abort("initialization from profile capability requires NetCDF");
#else
    m_init_uvtheta_profile = true;

    amrex::ParmParse pp_abl("ABL");
    std::string profileFile;
    pp_abl.query("init_profile", profileFile);

    auto ncf = ncutils::NCFile::open_par(
        profileFile, NC_NOWRITE | NC_NETCDF4 | NC_MPIIO,
        amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

    const auto num_prof_val = ncf.dim("nheight").len();

    m_theta_heights.resize(num_prof_val);
    m_theta_values.resize(num_prof_val);
    m_u_values.resize(num_prof_val);
    m_v_values.resize(num_prof_val);

    ncf.var("heights").get(m_theta_heights.data());
    ncf.var("theta").get(m_theta_values.data());
    ncf.var("u").get(m_u_values.data());
    ncf.var("v").get(m_v_values.data());

    m_thht_d.resize(num_prof_val);
    m_thvv_d.resize(num_prof_val);

    m_prof_u_d.resize(num_prof_val);
    m_prof_v_d.resize(num_prof_val);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_heights.begin(),
        m_theta_heights.end(), m_thht_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
        m_thvv_d.begin());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_u_values.begin(), m_u_values.end(),
        m_prof_u_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_v_values.begin(), m_v_values.end(),
        m_prof_v_d.begin());
#endif
}

void ABLFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& velocity,
    const amrex::Array4<amrex::Real>& density,
    const amrex::Array4<amrex::Real>& temperature) const
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    const amrex::Real rho_init = m_rho;

    const int ntvals = static_cast<int>(m_theta_heights.size());
    const int nwvals = static_cast<int>(m_wind_heights.size());
    const amrex::Real* th = m_thht_d.data();
    const amrex::Real* th_end = m_thht_d.end();
    const amrex::Real* tv = m_thvv_d.data();

    if (m_init_uvtheta_profile) {
        /*
         * Set wind and temperature profiles from netcdf input
         */
        const amrex::Real* uu = m_prof_u_d.data();
        const amrex::Real* vv = m_prof_v_d.data();

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            density(i, j, k) = rho_init;
            const auto idx = interp::bisection_search(th, th_end, z);
            temperature(i, j, k, 0) += interp::linear_impl(th, tv, z, idx);
            velocity(i, j, k, 0) += interp::linear_impl(th, uu, z, idx);
            velocity(i, j, k, 1) += interp::linear_impl(th, vv, z, idx);
        });
    } else if (m_initial_wind_profile) {
        const auto* xterrain_ptr = m_xterrain.data();
        const auto* yterrain_ptr = m_yterrain.data();
        const auto* zterrain_ptr = m_zterrain.data();
        const auto xterrain_size = m_xterrain.size();
        const auto yterrain_size = m_yterrain.size();
        const amrex::Real* windh = m_windht_d.data();
        const amrex::Real* uu = m_prof_u_d.data();
        const amrex::Real* vv = m_prof_v_d.data();
        const bool terrain_aligned_profile = m_terrain_aligned_profile;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real terrainHt =
                terrain_aligned_profile
                    ? interp::bilinear(
                          xterrain_ptr, xterrain_ptr + xterrain_size,
                          yterrain_ptr, yterrain_ptr + yterrain_size,
                          zterrain_ptr, x, y)
                    : 0.0_rt;
            z = amrex::max<amrex::Real>(0.5_rt * dx[2], z - terrainHt);
            density(i, j, k) = rho_init;
            const amrex::Real theta =
                (ntvals > 0) ? interp::linear(th, th + ntvals, tv, z) : tv[0];
            const amrex::Real umean_prof =
                (nwvals > 0) ? interp::linear(windh, windh + nwvals, uu, z)
                             : uu[0];
            const amrex::Real vmean_prof =
                (nwvals > 0) ? interp::linear(windh, windh + nwvals, vv, z)
                             : vv[0];

            temperature(i, j, k, 0) += theta;
            velocity(i, j, k, 0) += umean_prof;
            velocity(i, j, k, 1) += vmean_prof;
        });
    } else {
        /*
         * Set uniform/linear wind profile with specified temperature
         * profile
         */
        const bool linear_profile = m_linear_profile;

        const amrex::Real umean = !m_vel_timetable.empty()
                                      ? m_vel_speed * std::cos(m_vel_dir)
                                      : m_vel[0];
        const amrex::Real vmean = !m_vel_timetable.empty()
                                      ? m_vel_speed * std::sin(m_vel_dir)
                                      : m_vel[1];
        const amrex::Real wmean = !m_vel_timetable.empty() ? 0.0_rt : m_vel[2];

        const amrex::Real top_u_vel = m_top_vel[0];
        const amrex::Real top_v_vel = m_top_vel[1];
        const amrex::Real top_w_vel = m_top_vel[2];

        const amrex::Real bottom_u_vel = m_bottom_vel[0];
        const amrex::Real bottom_v_vel = m_bottom_vel[1];
        const amrex::Real bottom_w_vel = m_bottom_vel[2];

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            density(i, j, k) = rho_init;
            // Mean velocity field
            velocity(i, j, k, 0) = umean;
            velocity(i, j, k, 1) = vmean;
            velocity(i, j, k, 2) = wmean;

            amrex::Real theta = tv[0];
            for (int iz = 0; iz < ntvals - 1; ++iz) {
                if ((z > th[iz]) && (z <= th[iz + 1])) {
                    const amrex::Real slope =
                        (tv[iz + 1] - tv[iz]) / (th[iz + 1] - th[iz]);
                    theta = tv[iz] + ((z - th[iz]) * slope);
                }
            }

            temperature(i, j, k, 0) += theta;

            if (linear_profile) {
                velocity(i, j, k, 0) =
                    bottom_u_vel +
                    (z * (top_u_vel - bottom_u_vel) / (probhi[2] - problo[2]));
                velocity(i, j, k, 1) =
                    bottom_v_vel +
                    (z * (top_v_vel - bottom_v_vel) / (probhi[2] - problo[2]));
                velocity(i, j, k, 2) =
                    bottom_w_vel +
                    (z * (top_w_vel - bottom_w_vel) / (probhi[2] - problo[2]));
            }
        });
    }

    // velocity perturbations may be added on top of the simple wind
    // profiles specified in the input file or the general profiles from
    // a netcdf input
    if (m_perturb_vel) {
        const amrex::Real aval = m_Uperiods * 2.0_rt *
                                 std::numbers::pi_v<amrex::Real> /
                                 (probhi[1] - problo[1]);
        const amrex::Real bval = m_Vperiods * 2.0_rt *
                                 std::numbers::pi_v<amrex::Real> /
                                 (probhi[0] - problo[0]);
        const amrex::Real ufac = m_deltaU * std::exp(0.5_rt) / m_ref_height;
        const amrex::Real vfac = m_deltaV * std::exp(0.5_rt) / m_ref_height;
        const amrex::Real ref_height = m_ref_height;

        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real xl = x - problo[0];
            const amrex::Real yl = y - problo[1];
            const amrex::Real zl = z / ref_height;
            const amrex::Real damp = std::exp(-0.5_rt * zl * zl);

            velocity(i, j, k, 0) += ufac * damp * z * std::cos(aval * yl);
            velocity(i, j, k, 1) += vfac * damp * z * std::cos(bval * xl);
        });
    }
}

void ABLFieldInit::perturb_temperature(
    const int lev, const amrex::Geometry& geom, Field& temperature) const
{
    /** Perturbations for the temperature field
     *
     */

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    auto& theta_mf = temperature(lev);
    const auto theta_cutoff_height = m_theta_cutoff_height;
    const auto theta_gauss_mean = m_theta_gauss_mean;
    const auto theta_gauss_var = m_theta_gauss_var;
    const auto delta_t = m_delta_t;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(theta_mf, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& theta = theta_mf.array(mfi);

        amrex::ParallelForRNG(
            bx, [=] AMREX_GPU_DEVICE(
                    int i, int j, int k, const amrex::RandomEngine& engine) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                if (z < theta_cutoff_height) {
                    theta(i, j, k) = delta_t * amrex::RandomNormal(
                                                   theta_gauss_mean,
                                                   theta_gauss_var, engine);
                }
            });
    }
}

//! Initialize sfs tke field at the beginning of the simulation
void ABLFieldInit::init_tke(
    const amrex::Geometry& geom, amrex::MultiFab& tke_mf) const
{
    tke_mf.setVal(m_tke_init);

    if (!m_tke_init_profile && !m_initial_wind_profile) {
        return;
    }

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto tke_cutoff_height = m_tke_cutoff_height;
    const auto tke_init_factor = m_tke_init_factor;

    const auto& tke_arrs = tke_mf.arrays();
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    if (m_initial_wind_profile) {
        const auto* xterrain_ptr = m_xterrain.data();
        const auto* yterrain_ptr = m_yterrain.data();
        const auto* zterrain_ptr = m_zterrain.data();
        const amrex::Real* windh = m_windht_d.data();
        const bool terrain_aligned_profile = m_terrain_aligned_profile;
        const int nwvals = static_cast<int>(m_wind_heights.size());
        const amrex::Real* tke_data = m_prof_tke_d.data();
        const auto xterrain_size = m_xterrain.size();
        const auto yterrain_size = m_yterrain.size();
        amrex::ParallelFor(
            tke_mf, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                const amrex::Real terrainHt =
                    terrain_aligned_profile
                        ? interp::bilinear(
                              xterrain_ptr, xterrain_ptr + xterrain_size,
                              yterrain_ptr, yterrain_ptr + yterrain_size,
                              zterrain_ptr, x, y)
                        : 0.0_rt;
                z = amrex::max<amrex::Real>(0.5_rt * dx[2], z - terrainHt);
                const amrex::Real tke_prof =
                    (nwvals > 0)
                        ? interp::linear(windh, windh + nwvals, tke_data, z)
                        : tiny;

                tke_arrs[nbx](i, j, k) = tke_prof;
            });
    } else {
        // Profile definition from Beare et al. (2006)
        amrex::ParallelFor(
            tke_mf, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                if (z < tke_cutoff_height) {
                    tke_arrs[nbx](i, j, k) =
                        tke_init_factor *
                        std::pow(1.0_rt - (z / tke_cutoff_height), 3.0_rt);
                } else {
                    tke_arrs[nbx](i, j, k) = tiny;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
