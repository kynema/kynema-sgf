#include "amr-wind/utilities/tagging/BoxRefiner.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

// Adapted from OpenFOAM/src/meshTools/sets/cellSources/rotatedBoxToCell

namespace amr_wind::tagging {

namespace {

//! Utility function to parse inputs and return a vector instance
vs::Vector parse_vector(amrex::ParmParse& pp, const std::string& key)
{
    amrex::Vector<amrex::Real> tmp;
    pp.getarr(key, tmp);
    AMREX_ALWAYS_ASSERT(tmp.size() == 3U);

    return vs::Vector{tmp[0], tmp[1], tmp[2]};
}

/** Compute the corners that make up the hexahedral box
 *
 *  \return Position vectors of the 8 vertices that make up the box
 */
amrex::Vector<vs::Vector> compute_hex_corners(
    const vs::Vector& origin,
    const vs::Vector& x,
    const vs::Vector& y,
    const vs::Vector& z)
{

    amrex::Vector<vs::Vector> hex_nodes(8);
    hex_nodes[0] = origin;
    hex_nodes[1] = origin + x;
    hex_nodes[2] = origin + x + y;
    hex_nodes[3] = origin + y;
    hex_nodes[4] = origin + z;
    hex_nodes[5] = origin + z + x;
    hex_nodes[6] = origin + z + x + y;
    hex_nodes[7] = origin + z + y;

    return hex_nodes;
}

/** Compute the outward facing normals for the hexahedral box
 *
 *  \param hex_nodes Coordinates of the 8 vertices that make up the box
 *  \return Face normals for the six faces that make up the hexahedral box
 */
amrex::Vector<vs::Vector>
compute_face_normals(const amrex::Vector<vs::Vector>& hex_nodes)
{
    amrex::Vector<vs::Vector> face_normals(6);
    // Nodes that make up the face ordered such that the normal points outward
    amrex::Vector<amrex::Vector<int>> face_map{
        {0, 4, 7, 3}, // xmin
        {1, 2, 6, 5}, // xmax
        {0, 1, 5, 4}, // ymin
        {2, 3, 7, 6}, // ymax
        {0, 3, 2, 1}, // zmin
        {4, 5, 6, 7}  // zmax
    };

    for (int f = 0; f < 6; ++f) {
        const auto& fids = face_map[f];

        const vs::Vector v1 = hex_nodes[fids[1]] - hex_nodes[fids[0]];
        const vs::Vector v2 = hex_nodes[fids[3]] - hex_nodes[fids[0]];

        face_normals[f] = v1 ^ v2;
    }

    return face_normals;
}

} // namespace

BoxRefiner::BoxRefiner(const CFDSim& /*unused*/, const std::string& key)
{
    amrex::ParmParse pp(key);
    const auto origin = parse_vector(pp, "origin");
    const auto xaxis = parse_vector(pp, "xaxis");
    const auto yaxis = parse_vector(pp, "yaxis");
    const auto zaxis = parse_vector(pp, "zaxis");

    const auto hex_corners = compute_hex_corners(origin, xaxis, yaxis, zaxis);
    const auto face_normals = compute_face_normals(hex_corners);
    // Index of the starting corner for each face in the hex nodes list
    const amrex::Vector<int> face_origin{0, 1, 0, 2, 0, 4};

    // Define the bounding box around the hexahedral
    vs::Vector min_c(
        vs::DTraits<amrex::Real>::max(), vs::DTraits<amrex::Real>::max(),
        vs::DTraits<amrex::Real>::max());
    vs::Vector max_c(
        vs::DTraits<amrex::Real>::min(), vs::DTraits<amrex::Real>::min(),
        vs::DTraits<amrex::Real>::min());
    AMREX_ALWAYS_ASSERT(min_c.size() == hex_corners[0].size());
    AMREX_ALWAYS_ASSERT(max_c.size() == min_c.size());
    for (const auto& hc : hex_corners) {
        for (int i = 0; i < min_c.size(); i++) {
            min_c[i] = (min_c[i] > hc[i]) ? hc[i] : min_c[i];
            max_c[i] = (hc[i] > max_c[i]) ? hc[i] : max_c[i];
        }
    }
    vs::Vector search_radius = vs::Vector::zero();
    const amrex::Real search_fraction = 0.05_rt;
    AMREX_ALWAYS_ASSERT(search_radius.size() == min_c.size());
    for (int i = 0; i < min_c.size(); i++) {
        search_radius[i] = search_fraction * (max_c[i] - min_c[i]);
    }
    m_bound_box = amrex::RealBox(
        min_c[0] - search_radius[0], min_c[1] - search_radius[1],
        min_c[2] - search_radius[2], max_c[0] + search_radius[0],
        max_c[1] + search_radius[1], max_c[2] + search_radius[2]);

    // Setup data on device
    m_hex_corners.resize(8);
    m_face_normals.resize(6);
    m_face_origin.resize(6);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, hex_corners.begin(), hex_corners.end(),
        m_hex_corners.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, face_normals.begin(), face_normals.end(),
        m_face_normals.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, face_origin.begin(), face_origin.end(),
        m_face_origin.begin());
}

void BoxRefiner::operator()(
    const amrex::Box& bx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::TagBox::TagType>& tag) const
{
    const auto* hex_corners = m_hex_corners.data();
    const auto* face_normals = m_face_normals.data();
    const auto* fo = m_face_origin.data();
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
        const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
        const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

        // Position vector of cell center
        const vs::Vector pt(x, y, z);
        bool inside = true;

        // The point is inside the box if the signed distance is negative w.r.t.
        // all faces
        for (int f = 0; f < 6; ++f) {
            // Origin for this face
            const auto& forigin = hex_corners[fo[f]];
            // Cell center w.r.t. to face origin
            const auto ptloc = pt - forigin;

            if ((ptloc & face_normals[f]) > 0.0_rt) {
                inside = false;
            }
        }

        if (inside) {
            tag(i, j, k) = amrex::TagBox::SET;
        }
    });
}

} // namespace amr_wind::tagging
