#include "amr-wind/utilities/sampling/ProbeSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

ProbeSampler::ProbeSampler(const CFDSim& sim) : m_sim(sim) {}

ProbeSampler::~ProbeSampler() = default;

void ProbeSampler::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string pfile("probe_locations.txt");
    pp.query("probe_location_file", pfile);

    std::ifstream ifh(pfile, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find probe location file: " + pfile);
    }

    pp.queryarr("offsets", m_poffsets);
    if (!m_poffsets.empty()) {
        pp.getarr("offset_vector", m_offset_vector);
        AMREX_ALWAYS_ASSERT(
            static_cast<int>(m_offset_vector.size()) == AMREX_SPACEDIM);
    } else {
        // No offsets is implemented as 1 offset of 0.
        m_poffsets.push_back(0.0_rt);
    }

    int npts_file = 0;
    ifh >> npts_file;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    SampleLocType probes_file;
    // Read through points in file
    for (int i = 0; i < npts_file; ++i) {
        amrex::RealVect loc;
        ifh >> loc[0] >> loc[1] >> loc[2];
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        probes_file.push_back(loc, i);
    }
    AMREX_ALWAYS_ASSERT(probes_file.locations().size() == npts_file);
    // Incorporate offsets
    int idx = 0;
    m_npts = static_cast<int>(m_poffsets.size()) * npts_file;
    const auto& locs_file = probes_file.locations();
    for (amrex::Real m_poffset : m_poffsets) {
        for (int i = 0; i < npts_file; ++i) {
            amrex::RealVect loc;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                loc[d] = locs_file[i][d] + (m_poffset * m_offset_vector[d]);
            }
            m_probes.push_back(loc, idx);
            ++idx;
        }
    }
    AMREX_ALWAYS_ASSERT(m_probes.locations().size() == m_npts);

    check_bounds();
}

void ProbeSampler::check_bounds()
{
    const int lev = 0;
    const auto* prob_lo = m_sim.mesh().Geom(lev).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(lev).ProbHi();
    auto& probe_locs = m_probes.locations();
    bool all_ok = true;
    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (probe_locs[i][d] < (prob_lo[d] + bounds_tol)) {
                all_ok = false;
                probe_locs[i][d] = prob_lo[d] + (10 * bounds_tol);
            }
            if (probe_locs[i][d] > (prob_hi[d] - bounds_tol)) {
                all_ok = false;
                probe_locs[i][d] = prob_hi[d] - (10 * bounds_tol);
            }
        }
    }
    if (!all_ok) {
        amrex::Print() << "WARNING: ProbeSampler: Out of domain probe was "
                          "truncated to match domain"
                       << '\n';
    }
}

void ProbeSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const int lev = 0;
    const auto domain = m_sim.mesh().Geom(lev).Domain();
    sampling_locations(sample_locs, domain);
    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void ProbeSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const auto& probe_locs = m_probes.locations();
    const int lev = 0;
    const auto& dxinv = m_sim.mesh().Geom(lev).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(lev).ProbLoArray();
    for (int i = 0; i < m_npts; ++i) {
        const amrex::RealVect loc = {
            AMREX_D_DECL(probe_locs[i][0], probe_locs[i][1], probe_locs[i][2])};
        if (utils::contains(box, loc, plo, dxinv)) {
            sample_locs.push_back(loc, i);
        }
    }
}

#ifdef AMR_WIND_USE_NETCDF
void ProbeSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("offset_vector", m_offset_vector);
    grp.put_attr("offsets", m_poffsets);
}
#else
void ProbeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
#endif

} // namespace amr_wind::sampling
