#include "amr-wind/equation_systems/icns/source_terms/ABLMeanBoussinesq.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Boussinesq buoyancy source term for ABL simulations
 */
ABLMeanBoussinesq::ABLMeanBoussinesq(const CFDSim& sim)
    : m_mesh(sim.mesh()), m_transport(sim.transport_model())

{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_mean_boussinesq_term(this);

    // gravity in `incflo` namespace
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);

    // Backwards compatibility
    amrex::ParmParse pp_boussinesq_buoyancy("BoussinesqBuoyancy");
    amrex::ParmParse pp_abl("ABLMeanBoussinesq");

    bool read_temp_prof = false;
    if (pp_abl.contains("read_temperature_profile")) {
        pp_abl.get("read_temperature_profile", read_temp_prof);
        if (pp_boussinesq_buoyancy.contains("read_temperature_profile")) {
            amrex::Print()
                << "WARNING: BoussinesqBuoyancy.read_temperature_profile "
                   "option has been deprecated in favor of "
                   "ABLMeanBoussinesq.read_temperature_profile. Ignoring the"
                   "BoussinesqBuoyancy option in favor of the "
                   "ABLMeanBoussinesq "
                   "option."
                << '\n';
        }
    } else if (pp_boussinesq_buoyancy.contains("read_temperature_profile")) {
        amrex::Print()
            << "WARNING: BoussinesqBuoyancy.read_temperature_profile option "
               "has been deprecated in favor of "
               "ABLMeanBoussinesq.read_temperature_profile. Please replace "
               "this option."
            << '\n';
        pp_boussinesq_buoyancy.get("read_temperature_profile", read_temp_prof);
    }

    std::string tprofile_filename;
    if (pp_abl.contains("temperature_profile_filename")) {
        pp_abl.get("temperature_profile_filename", tprofile_filename);
        if (pp_boussinesq_buoyancy.contains("tprofile_filename")) {
            amrex::Print() << "WARNING: BoussinesqBuoyancy.tprofile_filename "
                              "option has been deprecated in favor of "
                              "ABLMeanBoussinesq.temperature_profile_filename. "
                              "Ignoring the"
                              "BoussinesqBuoyancy option in favor of the "
                              "ABLMeanBoussinesq "
                              "option."
                           << '\n';
        }
    } else if (pp_boussinesq_buoyancy.contains("tprofile_filename")) {
        amrex::Print()
            << "WARNING: BoussinesqBuoyancy.tprofile_filename option "
               "has been deprecated in favor of "
               "ABLMeanBoussinesq.temperature_profile_filename. Please replace "
               "this option."
            << '\n';
        pp_boussinesq_buoyancy.get("tprofile_filename", tprofile_filename);
    }

    if ((read_temp_prof) && (tprofile_filename.empty())) {
        m_const_profile = true;
        read_temperature_profile(tprofile_filename);
    } else {
        mean_temperature_init(abl.abl_statistics().theta_profile());
    }
}

ABLMeanBoussinesq::~ABLMeanBoussinesq() = default;

void ABLMeanBoussinesq::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    amrex::FArrayBox beta_fab(bx, 1, amrex::The_Async_Arena());
    amrex::Array4<amrex::Real> const& beta_arr = beta_fab.array();
    m_transport.beta_impl(lev, mfi, bx, beta_arr);

    amrex::FArrayBox ref_theta_fab(bx, 1, amrex::The_Async_Arena());
    amrex::Array4<amrex::Real> const& ref_theta_arr = ref_theta_fab.array();
    m_transport.ref_theta_impl(lev, mfi, bx, ref_theta_arr);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    // Mean temperature profile used to compute background forcing term

    const int idir = m_axis;
    const amrex::Real* theights = m_theta_ht.data();
    const amrex::Real* tvals = m_theta_vals.data();
    const amrex::Real* theights_end = m_theta_ht.end();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);
        const amrex::Real T0 = ref_theta_arr(i, j, k);
        const amrex::Real temp =
            amr_wind::interp::linear(theights, theights_end, tvals, ht);
        const amrex::Real fac = beta_arr(i, j, k) * (temp - T0);
        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

void ABLMeanBoussinesq::mean_temperature_init(const FieldPlaneAveraging& tavg)
{
    m_axis = tavg.axis();

    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(tavg.line_centroids().size()));
    m_theta_ht.resize(tavg.line_centroids().size());
    m_theta_vals.resize(tavg.line_average().size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_centroids().begin(),
        tavg.line_centroids().end(), m_theta_ht.begin());
    mean_temperature_update(tavg);
}

void ABLMeanBoussinesq::mean_temperature_update(const FieldPlaneAveraging& tavg)
{
    if (m_const_profile) {
        return;
    }
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tavg.line_average().begin(),
        tavg.line_average().end(), m_theta_vals.begin());
}

void ABLMeanBoussinesq::read_temperature_profile(
    const std::string& profile_file_name)
{

    amrex::Vector<amrex::Real> theta_ht, theta_vals;
    std::ifstream infile;
    int n_hts;
    infile.open(profile_file_name.c_str(), std::ios_base::in);
    infile >> n_hts;
    theta_ht.resize(n_hts);
    theta_vals.resize(n_hts);
    m_theta_ht.resize(n_hts);
    m_theta_vals.resize(n_hts);
    for (int i = 0; i < n_hts; i++) {
        infile >> theta_ht[i] >> theta_vals[i];
    }
    infile.close();

    // Now copy to GPU Device memory
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, theta_ht.begin(), theta_ht.end(),
        m_theta_ht.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, theta_vals.begin(), theta_vals.end(),
        m_theta_vals.begin());
}

} // namespace amr_wind::pde::icns
