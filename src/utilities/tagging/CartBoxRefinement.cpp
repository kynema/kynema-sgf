#include <iostream>
#include <fstream>

#include "amr-wind/utilities/tagging/CartBoxRefinement.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

namespace {

/** Read refinement bounding box boxes for a single level from a text file
 *
 *  The file must contain entries per line. The first entry is the total number
 *  of refinement boxes at the current level. This is followed by one entry per
 *  bbox that contains six floating point values. Corresponding to the
 *  coordinates of the min and max corners of the bounding box. The lines can
 *  contain text after the values and are ignored till EOF. For example,
 *
 *  ```
 *  4                               // Num. boxes in Lev 1
 *  -10.0 75.0 0.0 10.0 85.0 250.0  // T001
 *  -10.0 50.0 0.0 10.0 60.0 250.0  // T002
 *  -10.0 25.0 0.0 10.0 30.0 250.0  // T003
 *  -10.0 10.0 0.0 10.0 20.0 250.0  // T004
 *  ```
 *
 *  \param is Valid open input stream
 *  @return Vector containing RealBox instances for each bounding box
 */
amrex::Vector<amrex::RealBox> read_real_boxes(std::istream& is)
{
    amrex::Vector<amrex::RealBox> rbx_list;
    int nboxes;

    // Number of refinement regions at this level
    is >> nboxes;
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Read in the bounding boxes
    for (int ib = 0; ib < nboxes; ++ib) {
        amrex::Real xlo, ylo, zlo, xhi, yhi, zhi;
        is >> xlo >> ylo >> zlo >> xhi >> yhi >> zhi;
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        rbx_list.emplace_back(xlo, ylo, zlo, xhi, yhi, zhi);
    }

    return rbx_list;
}

amrex::BoxArray realbox_to_boxarray(
    const amrex::Vector<amrex::RealBox>& rbx, const amrex::Geometry& geom)
{
    amrex::BoxList bx_list;
    const auto* problo = geom.ProbLo();
    const auto* probhi = geom.ProbHi();
    const auto* dx = geom.CellSize();

    for (const auto& rb : rbx) {
        amrex::IntVect lo, hi;

        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            amrex::Real bbox_min = amrex::max(rb.lo()[i], problo[i]);
            amrex::Real bbox_max = amrex::min(rb.hi()[i], probhi[i]);
            amrex::Real rlo = std::floor((bbox_min - problo[i]) / dx[i]);
            amrex::Real rhi = std::ceil((bbox_max - problo[i]) / dx[i]);
            lo[i] = static_cast<int>(rlo);
            hi[i] = static_cast<int>(rhi);
        }
        bx_list.push_back({lo, hi});
    }

    return amrex::BoxArray(std::move(bx_list));
}

} // namespace

CartBoxRefinement::CartBoxRefinement(CFDSim& sim) : m_mesh(sim.mesh()) {}

void CartBoxRefinement::initialize(const std::string& key)
{
    std::string defn_file = "static_refinement.txt";
    {
        amrex::ParmParse pp(key);
        pp.query("static_refinement_def", defn_file);
    }

    std::ifstream ifh(defn_file, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find input file: " + defn_file);
    }

    read_inputs(m_mesh, ifh);
}

void CartBoxRefinement::read_inputs(
    const amrex::AmrCore& mesh, std::istream& ifh)
{
    const auto& geom = mesh.Geom();
    int max_lev = static_cast<int>(geom.size());

    int nlev_in;
    ifh >> nlev_in;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Issue a warning if the max levels in the input file is less than what's
    // requested in the refinement file.
    if (max_lev < nlev_in) {
        amrex::Print() << "WARNING: AmrMesh::finestLevel() is less than the "
                          "requested levels in static refinement file"
                       << '\n';
    }

    // Set the number of levels to the minimum of what is in the input file and
    // the simulation
    m_nlevels = amrex::min(max_lev, nlev_in);

    if (m_nlevels < 1) {
        return;
    }

    for (int lev = 0; lev < m_nlevels; ++lev) {
        auto rbx_list = read_real_boxes(ifh);
        auto ba = realbox_to_boxarray(rbx_list, geom[lev]);

        m_real_boxes.push_back(std::move(rbx_list));
        m_boxarrays.push_back(std::move(ba));
    }
}

void CartBoxRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    if (level < m_nlevels) {
        tags.setVal(m_boxarrays[level], amrex::TagBox::SET);
    }
}

} // namespace amr_wind
