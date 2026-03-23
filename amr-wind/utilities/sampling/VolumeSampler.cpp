#include <limits>

#include "amr-wind/utilities/sampling/VolumeSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::sampling {

VolumeSampler::VolumeSampler(const CFDSim& sim) : m_sim(sim) {}

VolumeSampler::~VolumeSampler() = default;

void VolumeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("hi", m_hi);
    pp.getarr("lo", m_lo);
    pp.getarr("num_points", m_npts_dir);
    check_bounds();
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_hi.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_lo.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == AMREX_SPACEDIM);

    // Update total number of points
    const size_t tmp =
        static_cast<size_t>(m_npts_dir[0]) * m_npts_dir[1] * m_npts_dir[2];
    if (tmp > static_cast<size_t>(std::numeric_limits<int>::max())) {
        amrex::Abort(
            "VolumeSampler: Volume definition " + key +
            " exceeds 32-bit integer limits");
    }
    m_npts = static_cast<int>(tmp);
}

void VolumeSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    bool all_ok = true;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_lo[d] < (prob_lo[d] + bounds_tol)) {
            all_ok = false;
            m_lo[d] = prob_lo[d] + (10 * bounds_tol);
        }
        if (m_lo[d] > (prob_hi[d] - bounds_tol)) {
            all_ok = false;
            m_lo[d] = prob_hi[d] - (10 * bounds_tol);
        }
        if (m_hi[d] < (prob_lo[d] + bounds_tol)) {
            all_ok = false;
            m_hi[d] = prob_lo[d] + (10 * bounds_tol);
        }
        if (m_hi[d] > (prob_hi[d] - bounds_tol)) {
            all_ok = false;
            m_hi[d] = prob_hi[d] - (10 * bounds_tol);
        }
    }
    if (!all_ok) {
        amrex::Print() << "WARNING: VolumeSampler: Out of domain corner was "
                          "truncated to match domain"
                       << '\n';
    }
}

void VolumeSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto domain = m_sim.mesh().Geom(lev).Domain();
    sampling_locations(sample_locs, domain);

    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void VolumeSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> dx = {
        (m_hi[0] - m_lo[0]) / m_npts_dir[0],
        (m_hi[1] - m_lo[1]) / m_npts_dir[1],
        (m_hi[2] - m_lo[2]) / m_npts_dir[2]};

    int idx = 0;
    for (int k = 0; k < m_npts_dir[2]; ++k) {
        for (int j = 0; j < m_npts_dir[1]; ++j) {
            for (int i = 0; i < m_npts_dir[0]; ++i) {
                const amrex::RealVect loc = {AMREX_D_DECL(
                    m_lo[0] + (dx[0] * i), m_lo[1] + (dx[1] * j),
                    m_lo[2] + (dx[2] * k))};
                if (utils::contains(box, loc, plo, dxinv)) {
                    sample_locs.push_back(loc, idx);
                }
                ++idx;
            }
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void VolumeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    const std::vector<int> ijk{m_npts_dir[0], m_npts_dir[1], m_npts_dir[2]};
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("ijk_dims", ijk);
    grp.put_attr("lo", m_lo);
    grp.put_attr("hi", m_hi);
}

void VolumeSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#else
void VolumeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void VolumeSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
