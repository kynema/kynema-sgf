#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLModulatedPowerLaw.H"
#include "amr-wind/wind_energy/ABLFillMPL.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_PlotFileUtil.H"
#include <sstream>
#include <iostream>
#include <string>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

ABLModulatedPowerLaw::ABLModulatedPowerLaw(CFDSim& sim)
    : m_sim(sim)
    , m_time(m_sim.time())
    , m_repo(m_sim.repo())
    , m_mesh(m_sim.mesh())
    , m_velocity(m_sim.repo().get_field("velocity"))
    , m_temperature(m_sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("MPL");

    pp.query("activate", m_activate_mpl);
    pp.query("zref", m_zref);
    pp.query("shear_exp", m_shear_exp);
    pp.query("umax_factor", m_umax_factor);
    pp.query("bulk_velocity", m_bulk_velocity);
    pp.query("shearlayer_height", m_shearlayer_height);
    pp.query("shearlayer_smear_thickness", m_shearlayer_smear_thickness);
    pp.query("wind_speed", m_wind_speed);
    pp.query("wind_direction", m_wind_direction);

    const amrex::Real wind_speed = m_wind_speed;
    const amrex::Real wind_direction = -m_wind_direction + 270.0_rt;
    const amrex::Real wind_direction_radian =
        amr_wind::utils::radians(wind_direction);

    m_uvec[0] = wind_speed * std::cos(wind_direction_radian);
    m_uvec[1] = wind_speed * std::sin(wind_direction_radian);
    m_uvec[2] = 0.0_rt;

    pp.query("start_time", m_start_time);
    pp.query("stop_time", m_stop_time);
    pp.query("degrees_per_second", m_degrees_per_sec);

    pp.query("delta_t", m_delta_t);
    pp.query("theta_cutoff_height", m_theta_cutoff_height);
    pp.query("theta_gauss_mean", m_theta_gauss_mean);
    pp.query("theta_gauss_var", m_theta_gauss_var);

    amrex::ParmParse pp_abl("ABL");
    pp_abl.getarr("temperature_heights", m_theta_heights);
    pp_abl.getarr("temperature_values", m_theta_values);

    AMREX_ALWAYS_ASSERT(m_theta_heights.size() == m_theta_values.size());
    int num_theta_values = static_cast<int>(m_theta_heights.size());
    m_thht_d.resize(num_theta_values);
    m_thvv_d.resize(num_theta_values);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_heights.begin(),
        m_theta_heights.end(), m_thht_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_theta_values.begin(), m_theta_values.end(),
        m_thvv_d.begin());
}

void ABLModulatedPowerLaw::post_init_actions()
{
    if (m_activate_mpl) {
        m_velocity.register_fill_patch_op<ABLFillMPL>(m_mesh, m_time, *this);
        m_temperature.register_fill_patch_op<ABLFillMPL>(m_mesh, m_time, *this);
    }
}

void ABLModulatedPowerLaw::pre_advance_work()
{

#ifdef AMR_WIND_USE_HELICS
    if (m_sim.helics().is_activated()) {
        const amrex::Real wind_speed =
            m_sim.helics().m_inflow_wind_speed_to_amrwind;
        const amrex::Real wind_direction =
            -m_sim.helics().m_inflow_wind_direction_to_amrwind + 270.0_rt;
        const amrex::Real wind_direction_radian =
            amr_wind::utils::radians(wind_direction);

        m_uvec[0] = wind_speed * std::cos(wind_direction_radian);
        m_uvec[1] = wind_speed * std::sin(wind_direction_radian);
        m_uvec[2] = 0.0_rt;
        return;
    }
#endif

    if (m_time.current_time() > m_start_time &&
        m_time.current_time() < m_stop_time) {
        m_wind_direction -= m_degrees_per_sec * m_time.delta_t();
    }
    const amrex::Real wind_direction = -m_wind_direction + 270.0_rt;
    const amrex::Real wind_direction_radian = utils::radians(wind_direction);
    m_uvec[0] = m_wind_speed * std::cos(wind_direction_radian);
    m_uvec[1] = m_wind_speed * std::sin(wind_direction_radian);
    m_uvec[2] = 0.0_rt;
}

void ABLModulatedPowerLaw::post_advance_work() {}

void ABLModulatedPowerLaw::set_velocity(
    const int lev,
    const amrex::Real /*time*/,
    const Field& fld,
    amrex::MultiFab& mfab,
    const int dcomp,
    const int orig_comp) const
{

    if (!m_activate_mpl) {
        return;
    }

    BL_PROFILE("amr-wind::ABLModulatedPowerLaw::set_velocity");

    const amrex::Real tvx = m_uvec[0];
    const amrex::Real tvy = m_uvec[1];
    const amrex::Real tvz = m_uvec[2];

    const amrex::Real zoffset = m_zoffset;
    const amrex::Real zref = m_zref;
    const amrex::Real shear_exp = m_shear_exp;
    const amrex::Real uref =
        vs::mag(vs::Vector{m_uvec[0], m_uvec[1], m_uvec[2]});
    const amrex::Real bulk_velocity = m_bulk_velocity;
    const amrex::Real umin = 0.0_rt;
    const amrex::Real umax_factor = m_umax_factor;
    const amrex::Real zc = m_shearlayer_height;
    const amrex::Real smear_coeff = 1.0_rt / m_shearlayer_smear_thickness;
    const amrex::Real z2 =
        (std::abs(shear_exp) > 0.0_rt)
            ? zref * std::pow(umax_factor, 1.0_rt / shear_exp)
            : 0.0_rt;
    const amrex::Real vmax = (z2 < zref) ? uref : uref * umax_factor;

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real height = probhi[2] - problo[2];
    const auto& dx = geom.CellSizeArray();

    const amrex::Real num1 = uref * std::pow(z2, shear_exp + 1.0_rt) /
                             (shear_exp + 1.0_rt) / std::pow(zref, shear_exp);
    const amrex::Real num2 = vmax * (height - z2);
    const amrex::Real denom =
        0.5_rt *
        (std::log(std::cosh(-smear_coeff * (height - zc))) / smear_coeff -
         std::log(std::cosh(-smear_coeff * (z2 - zc))) / smear_coeff +
         (height - z2));

    const amrex::Real upper_coeff =
        (bulk_velocity * height - num1 - num2) / denom;

    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow)) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            const auto& bx =
                utils::face_aware_boundary_box_intersection(gbx, dbx, ori);
            if (!bx.ok()) {
                continue;
            }

            const auto& arr = mfab[mfi].array();
            const int numcomp = mfab.nComp();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                const amrex::Real zeff = z - zoffset;
                amrex::Real pfac = (zeff > 0.0_rt)
                                       ? std::pow((zeff / zref), shear_exp)
                                       : 0.0_rt;
                pfac = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(pfac, umin), umax_factor);
                const amrex::Real tanhterm =
                    0.5_rt * upper_coeff *
                    (std::tanh(smear_coeff * (zeff - zc)) + 1.0_rt) / uref;

                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vels = {
                    AMREX_D_DECL(
                        tvx * (pfac + tanhterm), tvy * (pfac + tanhterm), tvz)};
                for (int n = 0; n < numcomp; n++) {
                    arr(i, j, k, dcomp + n) = vels[orig_comp + n];
                }
            });
        }
    }
}

void ABLModulatedPowerLaw::set_temperature(
    const int lev,
    const amrex::Real /*time*/,
    const Field& fld,
    amrex::MultiFab& mfab) const
{

    if (!m_activate_mpl) {
        return;
    }

    BL_PROFILE("amr-wind::ABLModulatedPowerLaw::set_temperature");

    const amrex::Real delta_t = m_delta_t;
    const amrex::Real theta_cutoff_height = m_theta_cutoff_height;
    const amrex::Real theta_gauss_mean = m_theta_gauss_mean;
    const amrex::Real theta_gauss_var = m_theta_gauss_var;

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    const int ntvals = static_cast<int>(m_theta_heights.size());
    const amrex::Real* th = m_thht_d.data();
    const amrex::Real* tv = m_thvv_d.data();

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow)) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            const auto& bx =
                utils::face_aware_boundary_box_intersection(gbx, dbx, ori);
            if (!bx.ok()) {
                continue;
            }

            const auto& arr = mfab[mfi].array();

            amrex::ParallelForRNG(
                bx,
                [=] AMREX_GPU_DEVICE(
                    int i, int j, int k, const amrex::RandomEngine& engine) {
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    amrex::Real theta = tv[0];
                    for (int iz = 0; iz < ntvals - 1; ++iz) {
                        if ((z > th[iz]) && (z <= th[iz + 1])) {
                            const amrex::Real slope =
                                (tv[iz + 1] - tv[iz]) / (th[iz + 1] - th[iz]);
                            theta = tv[iz] + ((z - th[iz]) * slope);
                        }
                    }
                    arr(i, j, k) = theta;

                    if (z < theta_cutoff_height) {
                        arr(i, j, k) += delta_t * amrex::RandomNormal(
                                                      theta_gauss_mean,
                                                      theta_gauss_var, engine);
                    }
                });
        }
    }
}

} // namespace amr_wind
