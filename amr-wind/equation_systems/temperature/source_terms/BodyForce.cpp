#include "amr-wind/equation_systems/temperature/source_terms/BodyForce.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_GpuContainers.H"
#include "AMReX_IntVect.H"
#include "AMReX_REAL.H"
#include "AMReX_Utility.H"
#include "AMReX_Vector.H"
#include <cstddef>
#include <fstream>
#include <ios>
#include <string>

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

BodyForce::BodyForce(const CFDSim& sim) : m_time(sim.time()), m_mesh(sim.mesh())
{

    amrex::ParmParse pp("BodyForceTemperature");
    pp.query("type", m_type);
    m_type = amrex::toLower(m_type);

    // Updated to underscores (more common) but still backwards-compatible
    if (m_type == "height_varying" || m_type == "height-varying") {
        pp.query("bodyforce-file", m_bforce_file);
        if (m_bforce_file.empty()) {
            pp.get("bodyforce_file", m_bforce_file);
        }
        read_bforce_profile(m_bforce_file);
    }
}
BodyForce::~BodyForce() = default;

void BodyForce::read_bforce_profile(const std::string& filename)
{
    std::ifstream infile;
    amrex::Vector<amrex::Real> bforce_theta;
    amrex::Vector<amrex::Real> bforce_hts;

    infile.open(filename.c_str(), std::ios_base::in);
    infile >> m_bforce_profile_nhts;

    bforce_theta.resize(m_bforce_profile_nhts);
    bforce_hts.resize(m_bforce_profile_nhts);

    m_prof_theta.resize(m_bforce_profile_nhts);
    m_ht.resize(m_bforce_profile_nhts);

    for (size_t i = 0; i < m_bforce_profile_nhts; i++) {
        infile >> bforce_hts[i] >> bforce_theta[i];
    }
    infile.close();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_hts.begin(), bforce_hts.end(),
        m_ht.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, bforce_theta.begin(), bforce_theta.end(),
        m_prof_theta.begin());
}

void BodyForce::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real* force_ht = m_ht.data();
    const amrex::Real* force_ht_end = m_ht.end();
    const amrex::Real* force_theta = m_prof_theta.data();

    if (m_type == "height_varying" || m_type == "height-varying") {

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::IntVect iv(i, j, k);
            const amrex::Real ht = problo[2] + ((iv[2] + 0.5_rt) * dx[2]);
            const amrex::Real ftheta = amr_wind::interp::linear(
                force_ht, force_ht_end, force_theta, ht);

            src_term(i, j, k, 0) += ftheta;
        });
    }
}

} // namespace amr_wind::pde::temperature
