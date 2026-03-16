#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <cstddef>
#include <ios>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Body Force
 */
BodyForce::BodyForce(const CFDSim& sim) : m_time(sim.time()), m_mesh(sim.mesh())
{

    // Body Force arguments
    amrex::ParmParse pp("BodyForce");
    pp.query("type", m_type);
    m_type = amrex::toLower(m_type);
    bool no_type_specified = !pp.contains("type");
    bool file_specified = pp.contains("uniform_timetable_file");

    // Prepare type of body force distribution
    if (m_type == "height_varying" || m_type == "height-varying") {
        // Constant in time, varies with z
        // Using underscores is preferred, remains backwards compatible
        pp.query("bodyforce-file", m_bforce_file);
        if (m_bforce_file.empty()) {
            pp.get("bodyforce_file", m_bforce_file);
        }
        read_bforce_profile(m_bforce_file);
    } else if (
        m_type == "uniform_timetable" ||
        (no_type_specified && file_specified)) {
        // Still used if type not specified but file is
        // Uniform in space, varies with time
        pp.get("uniform_timetable_file", m_utt_file);
        read_bforce_timetable(m_utt_file);
    } else {
        pp.getarr("magnitude", m_body_force);
        if (m_type == "oscillatory") {
            pp.get("angular_frequency", m_omega);
        } else if (m_type != "uniform_constant") {
            amrex::Abort(
                "BodyForce type not supported. Please choose uniform_constant "
                "(default), height_varying, oscillatory, or "
                "uniform_timetable.\n");
        }
    }
}

BodyForce::~BodyForce() = default;

void BodyForce::read_bforce_profile(const std::string& filename)
{

    std::ifstream infile;
    amrex::Vector<amrex::Real> bforce_x;
    amrex::Vector<amrex::Real> bforce_y;
    amrex::Vector<amrex::Real> bforce_hts;

    infile.open(filename.c_str(), std::ios_base::in);
    infile >> m_bforce_profile_nhts;

    bforce_x.resize(m_bforce_profile_nhts);
    bforce_y.resize(m_bforce_profile_nhts);
    bforce_hts.resize(m_bforce_profile_nhts);

    m_prof_x.resize(m_bforce_profile_nhts);
    m_prof_y.resize(m_bforce_profile_nhts);
    m_ht.resize(m_bforce_profile_nhts);
    for (size_t i = 0; i < m_bforce_profile_nhts; i++) {
        infile >> bforce_hts[i] >> bforce_x[i] >> bforce_y[i];
    }
    infile.close();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_x.begin(), bforce_x.end(),
        m_prof_x.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_y.begin(), bforce_y.end(),
        m_prof_y.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_hts.begin(), bforce_hts.end(),
        m_ht.begin());
}

void BodyForce::read_bforce_timetable(const std::string& filename)
{
    std::ifstream ifh(filename, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort(
            "Cannot find BodyForce uniform_timetable_file: " + filename + "\n");
    }
    amrex::Real data_time;
    amrex::Real data_fx;
    amrex::Real data_fy;
    amrex::Real data_fz;
    // Skip first line (header)
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_time) {
        ifh >> data_fx >> data_fy >> data_fz;
        m_time_table.push_back(data_time);
        m_fx_table.push_back(data_fx);
        m_fy_table.push_back(data_fy);
        m_fz_table.push_back(data_fz);
    }
}

void BodyForce::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& nph_time = 0.5_rt * (m_time.current_time() + m_time.new_time());
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real* force_ht = m_ht.data();
    const amrex::Real* force_ht_end = m_ht.end();
    const amrex::Real* force_x = m_prof_x.data();
    const amrex::Real* force_y = m_prof_y.data();

    if (m_type == "height_varying" || m_type == "height-varying") {

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::IntVect iv(i, j, k);
            const amrex::Real ht = problo[2] + ((iv[2] + 0.5_rt) * dx[2]);
            const amrex::Real fx =
                amr_wind::interp::linear(force_ht, force_ht_end, force_x, ht);
            const amrex::Real fy =
                amr_wind::interp::linear(force_ht, force_ht_end, force_y, ht);
            src_term(i, j, k, 0) += fx;
            src_term(i, j, k, 1) += fy;
        });

    } else {

        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
            m_body_force[0], m_body_force[1], m_body_force[2]};

        if (!m_utt_file.empty()) {
            // Populate forcing from file if supplied
            forcing[0] =
                amr_wind::interp::linear(m_time_table, m_fx_table, nph_time);
            forcing[1] =
                amr_wind::interp::linear(m_time_table, m_fy_table, nph_time);
            forcing[2] =
                amr_wind::interp::linear(m_time_table, m_fz_table, nph_time);
        }

        amrex::Real coeff =
            (m_type == "oscillatory") ? std::cos(m_omega * nph_time) : 1.0_rt;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            src_term(i, j, k, 0) += coeff * forcing[0];
            src_term(i, j, k, 1) += coeff * forcing[1];
            src_term(i, j, k, 2) += coeff * forcing[2];
        });
    }
}

} // namespace amr_wind::pde::icns
