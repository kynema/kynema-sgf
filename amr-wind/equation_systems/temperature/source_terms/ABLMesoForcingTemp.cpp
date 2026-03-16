#include "amr-wind/equation_systems/temperature/source_terms/ABLMesoForcingTemp.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Print.H"
#include <AMReX_REAL.H>
#include <memory>

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

ABLMesoForcingTemp::ABLMesoForcingTemp(const CFDSim& sim)
    : ABLMesoscaleForcing(sim, identifier())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_meso_temp_forcing(this);
    abl.abl_statistics().register_meso_temp_forcing(this);

    if (!abl.abl_meso_file().is_tendency_forcing()) {
        m_tendency = false;
        mean_temperature_init(
            abl.abl_statistics().theta_profile_fine(), abl.abl_meso_file());
    } else {
        m_tendency = true;
        mean_temperature_init(abl.abl_meso_file());
    }

    if ((amrex::toLower(m_forcing_scheme) == "indirect") &&
        !m_update_transition_height) {
        indirect_forcing_init(); // do this once
    }
}

ABLMesoForcingTemp::~ABLMesoForcingTemp() = default;

void ABLMesoForcingTemp::mean_temperature_init(const ABLMesoscaleInput& ncfile)
{
    const int num_meso_ht = ncfile.nheights();
    m_error_meso_avg_theta.resize(num_meso_ht);
    m_meso_theta_vals.resize(num_meso_ht);
    m_meso_ht.resize(num_meso_ht);
    m_err_Theta.resize(num_meso_ht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}
void ABLMesoForcingTemp::mean_temperature_init(
    const FieldPlaneAveragingFine& tavg, const ABLMesoscaleInput& ncfile)
{
    const int num_meso_ht = ncfile.nheights();
    m_nht = tavg.ncell_line();

    m_axis = tavg.axis();

    m_zht.resize(m_nht);
    m_theta_ht.resize(m_nht);
    m_err_Theta.resize(m_nht);
    m_error_meso_avg_theta.resize(m_nht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());

    std::copy(
        tavg.line_centroids().begin(), tavg.line_centroids().end(),
        m_zht.begin());

    m_meso_theta_vals.resize(num_meso_ht);
    m_meso_ht.resize(num_meso_ht);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, ncfile.meso_heights().begin(),
        ncfile.meso_heights().end(), m_meso_ht.begin());
}

amrex::Real ABLMesoForcingTemp::mean_temperature_heights(
    std::unique_ptr<ABLMesoscaleInput> const& ncfile)
{

    amrex::Real currtime;
    currtime = m_time.current_time();

    const amrex::Real interpTflux = amr_wind::interp::linear(
        ncfile->meso_times(), ncfile->meso_tflux(), currtime);

    if (m_forcing_scheme.empty()) {
        // no temperature profile assimilation
        return interpTflux;
    }

    const int num_meso_ht = ncfile->nheights();

    amrex::Vector<amrex::Real> time_interpolated_theta(num_meso_ht);

    for (int i = 0; i < num_meso_ht; i++) {
        time_interpolated_theta[i] = amr_wind::interp::linear(
            ncfile->meso_times(), ncfile->meso_temp(), currtime, num_meso_ht,
            i);
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, time_interpolated_theta.begin(),
        time_interpolated_theta.end(), m_error_meso_avg_theta.begin());

    for (int ih = 0; ih < num_meso_ht; ih++) {
        m_err_Theta[ih] = time_interpolated_theta[ih];
    }

    return interpTflux;
}

amrex::Real ABLMesoForcingTemp::mean_temperature_heights(
    const FieldPlaneAveragingFine& tavg,
    std::unique_ptr<ABLMesoscaleInput> const& ncfile)
{
    amrex::Real currtime;
    currtime = m_time.current_time();
    const auto& dt = m_time.delta_t();

    const amrex::Real interpTflux = amr_wind::interp::linear(
        ncfile->meso_times(), ncfile->meso_tflux(), currtime);

    if (m_forcing_scheme.empty()) {
        // no temperature profile assimilation
        return interpTflux;
    }

    amrex::Vector<amrex::Real> error_T(m_nht);

    for (int i = 0; i < m_nht; i++) {
        const amrex::Real interpolated_theta = amr_wind::interp::bilinear(

            ncfile->meso_times(), ncfile->meso_heights(), ncfile->meso_temp(),
            currtime, tavg.line_centroids()[i]);
        error_T[i] = interpolated_theta - tavg.line_average()[i];
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        if (m_update_transition_height) {
            // possible unexpected behaviors, as described in
            // ec5eb95c6ca853ce0fea8488e3f2515a2d6374e7
            // First the index in time

            m_transition_height = amr_wind::interp::linear(
                ncfile->meso_times(), ncfile->meso_transition_height(),
                currtime);
            amrex::Print() << "current transition height = "
                           << m_transition_height << '\n';

            set_transition_weighting();
            indirect_forcing_init();
        }

        amrex::Array<amrex::Real, 4> ezP_T;

        // form Z^T W y
        for (int i = 0; i < 4; i++) {
            ezP_T[i] = 0.0_rt;

            for (int ih = 0; ih < m_nht; ih++) {
                ezP_T[i] = ezP_T[i] + (error_T[ih] * m_W[ih] *
                                       std::pow(
                                           m_zht[ih] * m_scaleFact,
                                           static_cast<amrex::Real>(i)));
            }
        }

        for (int i = 0; i < 4; i++) {
            m_poly_coeff_theta[i] = 0.0_rt;
            for (int j = 0; j < 4; j++) {
                m_poly_coeff_theta[i] =
                    m_poly_coeff_theta[i] + (m_im_zTz(i, j) * ezP_T[j]);
            }
        }

        if (m_debug) {
            amrex::Print() << "direct vs indirect temperature error profile"
                           << '\n';
        }
        amrex::Vector<amrex::Real> error_T_direct(m_nht);
        for (int ih = 0; ih < m_nht; ih++) {
            error_T_direct[ih] = error_T[ih];
            error_T[ih] = 0.0_rt;
            for (int j = 0; j < 4; j++) {
                error_T[ih] = error_T[ih] + (m_poly_coeff_theta[j] *
                                             std::pow(
                                                 m_zht[ih] * m_scaleFact,
                                                 static_cast<amrex::Real>(j)));
            }

            if (m_debug) {
                amrex::Print() << m_zht[ih] << " " << error_T_direct[ih] << " "
                               << error_T[ih] << '\n';
            }
        }

        if (amrex::toLower(m_forcing_transition) == "indirecttodirect") {
            blend_forcings(error_T, error_T_direct, error_T);

            if (m_debug) {
                for (int ih = 0; ih < m_nht; ih++) {
                    amrex::Print() << m_zht[ih] << " " << error_T[ih] << '\n';
                }
            }
        }
    }

    if (forcing_to_constant()) {
        constant_forcing_transition(error_T);

        if (m_debug) {
            for (int ih = 0; ih < m_nht; ih++) {
                amrex::Print() << m_zht[ih] << " " << error_T[ih] << '\n';
            }
        }
    }

    for (int ih = 0; ih < m_nht; ih++) {
        error_T[ih] = error_T[ih] * m_gain_coeff / dt;
        m_err_Theta[ih] = error_T[ih];
    }

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, error_T.begin(), error_T.end(),
        m_error_meso_avg_theta.begin());

    return interpTflux;
}

void ABLMesoForcingTemp::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    if (m_forcing_scheme.empty()) {
        return;
    }

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    // The z values corresponding to the forcing values is either the number of
    // points in the netcdf input (for tendency) or the size of the plane
    // averaged temperature (non tendency)
    const amrex::Real* theights_begin =
        (m_tendency) ? m_meso_ht.data() : m_theta_ht.data();
    const amrex::Real* theights_end =
        (m_tendency) ? m_meso_ht.end() : m_theta_ht.end();
    const amrex::Real* theta_error_val = m_error_meso_avg_theta.data();

    const int idir = (int)m_axis;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);
        const amrex::Real theta_err = amr_wind::interp::linear(
            theights_begin, theights_end, theta_error_val, ht);

        // Compute Source term
        src_term(i, j, k, 0) += theta_err;
    });
}

} // namespace amr_wind::pde::temperature
