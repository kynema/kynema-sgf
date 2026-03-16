#include <limits>
#include "amr-wind/utilities/sampling/PlaneSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

PlaneSampler::PlaneSampler(const CFDSim& sim) : m_sim(sim) {}

PlaneSampler::~PlaneSampler() = default;

void PlaneSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("axis1", m_axis1);
    pp.getarr("axis2", m_axis2);
    pp.getarr("origin", m_origin);
    pp.getarr("num_points", m_npts_dir);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis1.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_axis2.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_npts_dir.size()) == 2);

    pp.queryarr("offsets", m_poffsets);
    int noffsets = static_cast<int>(m_poffsets.size());
    if (noffsets > 0) {
        if (pp.contains("normal")) {
            amrex::Abort(
                "PlaneSampler: option normal is deprecated and renamed to "
                "offset_vector");
        }

        pp.getarr("offset_vector", m_offset_vector);
        AMREX_ALWAYS_ASSERT(
            static_cast<int>(m_offset_vector.size()) == AMREX_SPACEDIM);
    } else {
        m_poffsets.push_back(0.0_rt);
    }

    check_bounds();

    // Update total number of points
    const size_t tmp = m_poffsets.size() * m_npts_dir[0] * m_npts_dir[1];
    if (tmp > static_cast<size_t>(std::numeric_limits<int>::max())) {
        amrex::Abort(
            "PlaneSampler: Plane definition " + key +
            " exceeds 32-bit integer limits");
    }
    m_npts = static_cast<int>(tmp);
}

void PlaneSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();

    // Do the checks based on a small fraction of the minimum dx
    const auto& fine_geom = m_sim.mesh().Geom(m_sim.mesh().finestLevel());
    amrex::Real min_dx = fine_geom.CellSize(0);
    for (int d = 1; d < AMREX_SPACEDIM; ++d) {
        min_dx = amrex::min<amrex::Real>(fine_geom.CellSize(d), min_dx);
    }
    const auto tol = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt * min_dx,
        bounds_tol);

    // First fix the origin so that it is within bounds, if it is close enough
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (amrex::Math::abs(m_origin[d] - prob_lo[d]) < tol) {
            m_origin[d] = prob_lo[d] + tol;
        }
        if (amrex::Math::abs(m_origin[d] - prob_hi[d]) < tol) {
            m_origin[d] = prob_hi[d] - tol;
        }
    }

    // Fix the axis so that it is within bounds, if it is close enough
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (amrex::Math::abs(m_origin[d] + m_axis1[d] - prob_lo[d]) <
            2.0_rt * tol) {
            m_axis1[d] = prob_lo[d] - m_origin[d] + (2.0_rt * tol);
        }
        if (amrex::Math::abs(m_origin[d] + m_axis1[d] - prob_hi[d]) <
            2.0_rt * tol) {
            m_axis1[d] = prob_hi[d] - m_origin[d] - (2.0_rt * tol);
        }
        if (amrex::Math::abs(m_origin[d] + m_axis2[d] - prob_lo[d]) <
            2.0_rt * tol) {
            m_axis2[d] = prob_lo[d] - m_origin[d] + (2.0_rt * tol);
        }
        if (amrex::Math::abs(m_origin[d] + m_axis2[d] - prob_hi[d]) <
            2.0_rt * tol) {
            m_axis2[d] = prob_hi[d] - m_origin[d] - (2.0_rt * tol);
        }
    }

    const int nplanes = static_cast<int>(m_poffsets.size());
    for (int k = 0; k < nplanes; ++k) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const auto point =
                m_origin[d] + (m_poffsets[k] * m_offset_vector[d]);
            const amrex::Vector<amrex::Real> points = {
                point, point + m_axis1[d], point + m_axis2[d]};
            for (const auto& pt : points) {
                if ((pt < prob_lo[d]) || (pt >= prob_hi[d])) {
                    amrex::Abort(
                        "PlaneSampler: Point out of domain. Redefine your "
                        "planes so they are completely inside the domain.");
                }
            }
        }
    }
}

void PlaneSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto domain = m_sim.mesh().Geom(lev).Domain();
    sampling_locations(sample_locs, domain);

    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void PlaneSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    amrex::Array<amrex::Real, AMREX_SPACEDIM> dx;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> dy;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dx[d] = m_axis1[d] / amrex::max(m_npts_dir[0] - 1, 1);
        dy[d] = m_axis2[d] / amrex::max(m_npts_dir[1] - 1, 1);
    }

    int idx = 0;
    const int nplanes = static_cast<int>(m_poffsets.size());
    const int lev = 0;
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    for (int k = 0; k < nplanes; ++k) {
        for (int j = 0; j < m_npts_dir[1]; ++j) {
            for (int i = 0; i < m_npts_dir[0]; ++i) {
                amrex::RealVect loc;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    loc[d] = m_origin[d] + (dx[d] * i) + (dy[d] * j) +
                             (m_poffsets[k] * m_offset_vector[d]);
                }
                if (utils::contains(box, loc, plo, dxinv)) {
                    sample_locs.push_back(loc, idx);
                }
                ++idx;
            }
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void PlaneSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    const std::vector<int> ijk{
        m_npts_dir[0], m_npts_dir[1], static_cast<int>(m_poffsets.size())};
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("ijk_dims", ijk);
    grp.put_attr("origin", m_origin);
    grp.put_attr("axis1", m_axis1);
    grp.put_attr("axis2", m_axis2);
    grp.put_attr("offset_vector", m_offset_vector);
    grp.put_attr("offsets", m_poffsets);
}

void PlaneSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#else
void PlaneSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void PlaneSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
