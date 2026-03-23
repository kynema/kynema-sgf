#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/index_operations.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

LineSampler::LineSampler(const CFDSim& sim) : m_sim(sim) {}

LineSampler::~LineSampler() = default;

void LineSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.get("num_points", m_npts);
    pp.getarr("start", m_start);
    pp.getarr("end", m_end);

    check_bounds();
}

void LineSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    bool all_ok = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_start[d] < (prob_lo[d] + bounds_tol)) {
            all_ok = false;
            m_start[d] = prob_lo[d] + (10 * bounds_tol);
        }
        if (m_start[d] > (prob_hi[d] - bounds_tol)) {
            all_ok = false;
            m_start[d] = prob_hi[d] - (10 * bounds_tol);
        }
        if (m_end[d] < (prob_lo[d] + bounds_tol)) {
            all_ok = false;
            m_end[d] = prob_lo[d] + (10 * bounds_tol);
        }
        if (m_end[d] > (prob_hi[d] - bounds_tol)) {
            all_ok = false;
            m_end[d] = prob_hi[d] - (10 * bounds_tol);
        }
    }
    if (!all_ok) {
        amrex::Print() << "WARNING: LineSampler: Out of domain line was "
                          "truncated to match domain"
                       << '\n';
    }
}

void LineSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto domain = m_sim.mesh().Geom(lev).Domain();
    sampling_locations(sample_locs, domain);

    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void LineSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    const amrex::Real ndiv = amrex::max(m_npts - 1, 1);
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = (m_end[d] - m_start[d]) / ndiv;
    }

    for (int i = 0; i < m_npts; ++i) {
        const amrex::RealVect loc = {AMREX_D_DECL(
            m_start[0] + (i * dx[0]), m_start[1] + (i * dx[1]),
            m_start[2] + (i * dx[2]))};
        if (utils::contains(box, loc, plo, dxinv)) {
            sample_locs.push_back(loc, i);
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void LineSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
}

void LineSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#else
void LineSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void LineSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
