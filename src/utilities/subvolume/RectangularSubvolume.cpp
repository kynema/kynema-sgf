#include "src/utilities/subvolume/RectangularSubvolume.H"
#include "src/utilities/constants.H"
#include "src/CFDSim.H"

#include <cmath>

#include "AMReX_ParmParse.H"

namespace kynema_sgf::subvolume {

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
    constexpr amrex::Real tol = 1.0e-4_rt;

    if (m_npts_vec.size() != AMREX_SPACEDIM ||
        m_dx_vec.size() != AMREX_SPACEDIM ||
        m_origin.size() != AMREX_SPACEDIM) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": origin, num_points, and dx_vec/dx must all have " +
            std::to_string(AMREX_SPACEDIM) + " entries.");
    }

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (m_npts_vec[d] <= 0) {
            amrex::Abort(
                "RectangularSubvolume " + m_label +
                ": num_points entries must all be positive.");
        }
        if (m_dx_vec[d] <= 0.0_rt) {
            amrex::Abort(
                "RectangularSubvolume " + m_label +
                ": dx and dx_vec entries must all be positive.");
        }
    }

    bool found = false;
    const auto& geom = m_sim.mesh().Geom();
    amrex::IntVect best_start(0, 0, 0);
    amrex::IntVect best_stride(1, 1, 1);
    amrex::Box best_src_box;

    for (int i = 0; i < m_sim.repo().num_active_levels(); i++) {
        amrex::IntVect stride(1, 1, 1);
        bool valid_stride = true;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const auto ratio = m_dx_vec[d] / geom[i].CellSize(d);
            const int iratio = static_cast<int>(std::lround(ratio));
            if (iratio < 1 ||
                std::abs(ratio - static_cast<amrex::Real>(iratio)) > tol) {
                valid_stride = false;
                break;
            }
            stride[d] = iratio;
        }

        if (!valid_stride) {
            continue;
        }

        int i0 = static_cast<int>(std::lround(
            (m_origin[0] - geom[i].ProbLo(0)) / geom[i].CellSize(0)));
        int j0 = static_cast<int>(std::lround(
            (m_origin[1] - geom[i].ProbLo(1)) / geom[i].CellSize(1)));
        int k0 = static_cast<int>(std::lround(
            (m_origin[2] - geom[i].ProbLo(2)) / geom[i].CellSize(2)));

        const bool aligned_origin =
            std::abs(
                geom[i].ProbLo(0) + i0 * geom[i].CellSize(0) - m_origin[0]) <
                tol &&
            std::abs(
                geom[i].ProbLo(1) + j0 * geom[i].CellSize(1) - m_origin[1]) <
                tol &&
            std::abs(
                geom[i].ProbLo(2) + k0 * geom[i].CellSize(2) - m_origin[2]) <
                tol;

        if (!aligned_origin) {
            continue;
        }

        amrex::Box bx(
            amrex::IntVect(i0, j0, k0),
            amrex::IntVect(
                i0 + (m_npts_vec[0] - 1) * stride[0],
                j0 + (m_npts_vec[1] - 1) * stride[1],
                k0 + (m_npts_vec[2] - 1) * stride[2]));

        if (!m_sim.mesh().boxArray()[i].contains(bx)) {
            continue;
        }

        found = true;
        m_lev_for_sub = i;
        best_start = amrex::IntVect(i0, j0, k0);
        best_stride = stride;
        best_src_box = bx;
    }

    if (!found) {
        amrex::Abort(
            "RectangularSubvolume " + m_label +
            ": Could not map requested output mesh to an existing level. "
            "The requested dx or dx_vec must be an integer multiple of a level "
            "cell size in each direction, the origin must align to that "
            "level grid, and the full sampled extent must be contained in "
            "that level.");
    }

    m_stride = best_stride;
    if (m_verbose > 0) {
        amrex::Print() << "RectangularSubvolume " + m_label
                       << ": Using source level " << m_lev_for_sub
                       << " with output stride (" << m_stride[0] << ", "
                       << m_stride[1] << ", " << m_stride[2] << ")\n"
                       << "RectangularSubvolume " + m_label
                       << ": Source sampling extent is " << best_src_box
                       << "\n";
    }

    amrex::Box out_box(
        amrex::IntVect(0, 0, 0),
        amrex::IntVect(
            m_npts_vec[0] - 1, m_npts_vec[1] - 1, m_npts_vec[2] - 1));

    if (!m_chunk_size_present && m_chunk_size_vec.empty()) {
        m_chunk_size_vec.resize(AMREX_SPACEDIM);
        m_chunk_size_vec[0] = m_sim.mesh().maxGridSize(m_lev_for_sub)[0];
        m_chunk_size_vec[1] = m_sim.mesh().maxGridSize(m_lev_for_sub)[1];
        m_chunk_size_vec[2] = m_sim.mesh().maxGridSize(m_lev_for_sub)[2];
    }

    amrex::IntVect chunk_size(
        m_chunk_size_vec[0], m_chunk_size_vec[1], m_chunk_size_vec[2]);

    amrex::BoxArray ba(out_box);
    ba.maxSize(chunk_size);

    // Create a new box array: same layout but different box definitions
    amrex::BoxArray ba_src(ba); // shallow copy to get same layout, size
    ba_src.uniqify();           // force deep copy to avoid changing original
    for (int ib = 0; ib < ba_src.size(); ++ib) {
        const auto b = ba_src[ib];
        const amrex::IntVect lo(
            best_start[0] + b.smallEnd(0) * m_stride[0],
            best_start[1] + b.smallEnd(1) * m_stride[1],
            best_start[2] + b.smallEnd(2) * m_stride[2]);
        const amrex::IntVect hi(
            best_start[0] + b.bigEnd(0) * m_stride[0],
            best_start[1] + b.bigEnd(1) * m_stride[1],
            best_start[2] + b.bigEnd(2) * m_stride[2]);
        ba_src.set(ib, amrex::Box(lo, hi));
    }

    const amrex::Real dx0 = m_stride[0] * geom[m_lev_for_sub].CellSize(0);
    const amrex::Real dx1 = m_stride[1] * geom[m_lev_for_sub].CellSize(1);
    const amrex::Real dx2 = m_stride[2] * geom[m_lev_for_sub].CellSize(2);
    // Normalize to the exact sampled spacing so output geometry remains
    // consistent with the integer stride mapping.
    m_dx_vec[0] = dx0;
    m_dx_vec[1] = dx1;
    m_dx_vec[2] = dx2;

    amrex::RealBox out_real_box(
        {m_origin[0] + 0.5_rt * (geom[m_lev_for_sub].CellSize(0) - dx0),
         m_origin[1] + 0.5_rt * (geom[m_lev_for_sub].CellSize(1) - dx1),
         m_origin[2] + 0.5_rt * (geom[m_lev_for_sub].CellSize(2) - dx2)},
        {m_origin[0] + 0.5_rt * (geom[m_lev_for_sub].CellSize(0) - dx0) +
             m_npts_vec[0] * dx0,
         m_origin[1] + 0.5_rt * (geom[m_lev_for_sub].CellSize(1) - dx1) +
             m_npts_vec[1] * dx1,
         m_origin[2] + 0.5_rt * (geom[m_lev_for_sub].CellSize(2) - dx2) +
             m_npts_vec[2] * dx2});

    m_out_geom = amrex::Geometry(
        out_box, out_real_box, geom[m_lev_for_sub].Coord(),
        geom[m_lev_for_sub].isPeriodic());

    if (m_verbose > 0) {
        amrex::Print() << "RectangularSubvolume " + m_label + ": BoxArray is "
                       << ba << "\n";
    }

    m_ba = ba;
    m_src_ba = ba_src;
}

} // namespace kynema_sgf::subvolume
