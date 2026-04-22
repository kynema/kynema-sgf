#include "src/utilities/subvolume/RectangularSubvolume.H"
#include "src/utilities/constants.H"
#include "src/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::subvolume {

RectangularSubvolume::RectangularSubvolume(const CFDSim& sim) : m_sim(sim) {}

RectangularSubvolume::~RectangularSubvolume() = default;

void RectangularSubvolume::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    pp.getarr("origin", m_origin);
    pp.getarr("num_points", m_npts_vec);
    amrex::Real dx{1.0_rt};
    pp.query("dx", dx);
    if (pp.contains("dx")) {
        m_dx_vec.resize(AMREX_SPACEDIM);
        m_dx_vec[0] = dx;
        m_dx_vec[1] = dx;
        m_dx_vec[2] = dx;
    } else {
        pp.getarr("dx_vec", m_dx_vec);
    }

    pp.queryarr("chunk_size_vec", m_chunk_size_vec);
    m_chunk_size_present = pp.contains("chunk_size_vec");

    pp.query("verbose", m_verbose);
}

void RectangularSubvolume::evaluate_inputs()
{
    bool found = false;
    const auto& geom = m_sim.mesh().Geom();
    for (int i = 0; i < m_sim.repo().num_active_levels(); i++) {
        if (!found) {
            if (std::abs(m_dx_vec[0] - geom[i].CellSize(0)) < 1e-4 &&
                std::abs(m_dx_vec[1] - geom[i].CellSize(1)) < 1e-4 &&
                std::abs(m_dx_vec[2] - geom[i].CellSize(2)) < 1e-4) {

                amrex::Print()
                    << "RectangularSubvolume " + m_label +
                           ": Resolution specified matches that of level "
                    << i << "\n";
                found = true;
                m_lev_for_sub = i;
            }
        }
    }

    if (!found) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Resolution specified for subvolume does not match the "
            "resolution of any of the mesh levels.");
    }

    // **************************************************************
    // Now that we know which level we're at, we can figure out which
    // (i,j,k) the origin corresponds to
    // **************************************************************
    int i0 = static_cast<int>(std::lround(
        (m_origin[0] - geom[m_lev_for_sub].ProbLo(0)) /
        geom[m_lev_for_sub].CellSize(0)));
    int j0 = static_cast<int>(std::lround(
        (m_origin[1] - geom[m_lev_for_sub].ProbLo(1)) /
        geom[m_lev_for_sub].CellSize(1)));
    int k0 = static_cast<int>(std::lround(
        (m_origin[2] - geom[m_lev_for_sub].ProbLo(2)) /
        geom[m_lev_for_sub].CellSize(2)));

    found = false;
    if (std::abs(
            geom[m_lev_for_sub].ProbLo(0) +
            i0 * geom[m_lev_for_sub].CellSize(0) - m_origin[0]) < 1e-4 &&
        std::abs(
            geom[m_lev_for_sub].ProbLo(1) +
            j0 * geom[m_lev_for_sub].CellSize(1) - m_origin[1]) < 1e-4 &&
        std::abs(
            geom[m_lev_for_sub].ProbLo(2) +
            k0 * geom[m_lev_for_sub].CellSize(2) - m_origin[2]) < 1e-4) {
        amrex::Print()
            << "RectangularSubvolume " + m_label +
                   ": Specified origin is the lower left corner of cell "
            << amrex::IntVect(i0, j0, k0) << "\n";
        found = true;
    }

    if (!found) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Origin specified does not correspond to a node at this level.");
    }

    amrex::Box bx(
        amrex::IntVect(i0, j0, k0),
        amrex::IntVect(
            i0 + m_npts_vec[0] - 1, j0 + m_npts_vec[1] - 1,
            k0 + m_npts_vec[2] - 1));
    if (m_verbose > 0) {
        amrex::Print() << "RectangularSubvolume " + m_label +
                              ": Box requested is "
                       << bx << "\n";
    }

    if (!m_sim.mesh().boxArray()[m_lev_for_sub].contains(bx)) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Box requested is larger than the existing domain");
    }

    if (!m_chunk_size_present && m_chunk_size_vec.empty()) {
        m_chunk_size_vec.resize(AMREX_SPACEDIM);
        m_chunk_size_vec[0] = m_sim.mesh().maxGridSize(m_lev_for_sub)[0];
        m_chunk_size_vec[1] = m_sim.mesh().maxGridSize(m_lev_for_sub)[1];
        m_chunk_size_vec[2] = m_sim.mesh().maxGridSize(m_lev_for_sub)[2];
    }

    amrex::IntVect chunk_size(
        m_chunk_size_vec[0], m_chunk_size_vec[1], m_chunk_size_vec[2]);

    amrex::BoxArray ba(bx);
    ba.maxSize(chunk_size);

    if (m_verbose > 0) {
        amrex::Print() << "RectangularSubvolume " + m_label + ": BoxArray is "
                       << ba << "\n";
    }

    m_ba = ba;
}

} // namespace amr_wind::subvolume
