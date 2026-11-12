#ifndef TEST_EXT_TURBINE_UTILS_H
#define TEST_EXT_TURBINE_UTILS_H

#include <iterator>
#include "gtest/gtest.h"
#include "AMReX_ParallelDescriptor.H"

// The external turbine data hold an MPI communicator; without MPI, AMReX
// provides stand-ins
#ifndef BL_USE_MPI
using amrex::MPI_Comm;
using amrex::MPI_COMM_NULL;
#endif

#include "src/actuator/turbine/external/turbine_external_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

/** Grid and metadata of an external turbine, without the external solver
 *
 *  Holds what make_component_views needs to build the hub, blade and tower
 *  views.
 */
template <typename SolverTurbine, typename SolverData>
struct ExtTurbineViewData
{
    using GridType = kynema_sgf::actuator::ActGrid;
    using MetaType =
        kynema_sgf::actuator::TurbineExternalData<SolverTurbine, SolverData>;

    GridType& grid() { return m_grid; }
    MetaType& meta() { return m_meta; }

    GridType m_grid;
    MetaType m_meta;
};

/** Fill the point lists of a turbine in the order OpenFAST gives them
 *
 *  Force points: the hub, the nodes of each blade from root to tip, the
 *  tower. Velocity points: the same order, hub first, one for each force
 *  point and at the same position. The position tells which point it is: the
 *  hub is at the origin, node k of blade ib is at (0, ib + 1, k + 1) and node
 *  j of the tower is at (0, 0, -(j + 1)).
 */
template <typename datatype>
void init_hub_first_points(
    datatype& data,
    const int num_blades,
    const int num_pts_blade,
    const int num_pts_tower)
{
    namespace vs = kynema_sgf::vs;
    auto& grid = data.grid();
    auto& meta = data.meta();
    meta.num_blades = num_blades;
    meta.num_pts_blade = num_pts_blade;
    meta.num_vel_pts_blade = num_pts_blade;
    meta.num_pts_tower = num_pts_tower;

    const int npts = 1 + num_blades * num_pts_blade + num_pts_tower;
    grid.resize(npts, npts);
    meta.chord.assign(npts, 1.0_rt);
    meta.vel_rel.assign(npts, vs::Vector::zero());

    grid.pos[0] = vs::Vector::zero();
    for (int ib = 0; ib < num_blades; ++ib) {
        for (int k = 0; k < num_pts_blade; ++k) {
            const int ip = 1 + ib * num_pts_blade + k;
            grid.pos[ip] = vs::Vector{
                0.0_rt, static_cast<amrex::Real>(ib + 1),
                static_cast<amrex::Real>(k + 1)};
        }
    }
    for (int j = 0; j < num_pts_tower; ++j) {
        const int ip = 1 + num_blades * num_pts_blade + j;
        grid.pos[ip] =
            vs::Vector{0.0_rt, 0.0_rt, -static_cast<amrex::Real>(j + 1)};
    }
    for (int ip = 0; ip < npts; ++ip) {
        grid.vel_pos[ip] = grid.pos[ip];
    }
}

/** Check that every blade view pairs force point k with its velocity point
 *
 *  Entry k of the velocity views of a blade must be the velocity point of
 *  force point k of that blade. The filtered lifting line correction relies
 *  on it: it adds the correction computed for force point k to velocity
 *  entry k.
 */
template <typename datatype>
void check_blade_views_aligned(datatype& data)
{
    namespace vs = kynema_sgf::vs;
    const auto& grid = data.grid();
    const auto& meta = data.meta();
    const int npb = meta.num_pts_blade;

    ASSERT_EQ(static_cast<int>(meta.blades.size()), meta.num_blades);
    for (int ib = 0; ib < meta.num_blades; ++ib) {
        const auto& cv = meta.blades[ib];
        const int start = 1 + ib * npb;
        EXPECT_EQ(std::distance(grid.pos.data(), cv.pos.data()), start);
        EXPECT_EQ(std::distance(grid.vel.data(), cv.vel.data()), start);
        EXPECT_EQ(
            std::distance(grid.vel_pos.data(), cv.vel_pos.data()), start);
        ASSERT_EQ(static_cast<int>(cv.vel_pos.size()), npb);
        for (int k = 0; k < npb; ++k) {
            EXPECT_EQ(vs::mag(cv.vel_pos[k] - cv.pos[k]), 0.0_rt)
                << "blade " << ib << ", node " << k;
        }
    }
}

} // namespace kynema_sgf_tests

#endif /* TEST_EXT_TURBINE_UTILS_H */
