#include "ks_test_utils/MeshTest.H"
#include "test_ext_turbine_utils.H"

#include "src/actuator/ActuatorContainer.H"
#include "src/actuator/turbine/external/turbine_external_ops.H"
#include "src/core/vs/vector_space.H"

#include <limits>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

struct MockBoundsTurbine
{};

struct MockBoundsSolverData
{};

using MockBoundsData =
    ExtTurbineViewData<MockBoundsTurbine, MockBoundsSolverData>;

//! Container with access to the positions and velocities of its points
class BoundsActContainer : public kynema_sgf::actuator::ActuatorContainer
{
public:
    BoundsActContainer(amrex::AmrCore& mesh, const int nobj)
        : kynema_sgf::actuator::ActuatorContainer(mesh, nobj)
    {}

    kynema_sgf::actuator::ActuatorCloud& get_data_obj() { return point_data(); }
};

//! Domain of 128 in every direction, with walls in z
class ActuatorBoundsTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{128.0_rt, 128.0_rt, 128.0_rt}};
            amrex::Vector<int> periodic{{1, 1, 0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", periodic);
        }
    }
};

} // namespace

TEST(ActuatorBounds, gaussian_beyond_bound_box_is_counted)
{
    namespace vs = kynema_sgf::vs;
    namespace external = kynema_sgf::actuator::external;
    const int num_blades = 3;
    const int npb = 5;
    const int npt = 4;
    MockBoundsData data;
    init_hub_first_points(data, num_blades, npb, npt);
    auto& meta = data.meta();
    external::make_component_views(data);

    // Box made from the inputs of a turbine at the origin with a rotor
    // radius of 8. The points are up to 5 away from the hub. init_epsilon
    // adds the reach of the Gaussians, 4 epsilon = 2, to the box
    data.info().bound_box =
        amrex::RealBox(-10.0_rt, -10.0_rt, -10.0_rt, 10.0_rt, 10.0_rt, 10.0_rt);
    meta.eps_inp = vs::Vector{0.5_rt, 0.5_rt, 0.5_rt};
    external::init_epsilon(data);
    EXPECT_EQ(external::check_bound_box(data), 0);
    EXPECT_FALSE(meta.bound_box_warned);

    // The turbine moves by 11 along x: the Gaussians of the points of the
    // blades and of the tower now reach x = 13, beyond the face of the box
    // at 12. The hub has no Gaussian (no nacelle drag in this turbine)
    for (auto& pos : data.grid().pos) {
        pos.x() += 11.0_rt;
    }
    EXPECT_EQ(external::check_bound_box(data), num_blades * npb + npt);
    EXPECT_TRUE(meta.bound_box_warned);
}

TEST_F(ActuatorBoundsTest, points_outside_domain_are_counted)
{
    namespace vs = kynema_sgf::vs;
    initialize_mesh();
    auto& vel = sim().repo().declare_field("velocity", 3, 3);
    auto& density = sim().repo().declare_field("density", 1, 3);
    vel.setVal(2.0_rt);
    density.setVal(1.0_rt);

    const int num_pts = 6;
    BoundsActContainer ac(mesh(), 1);
    auto& data = ac.get_data_obj();
    data.num_pts[0] = num_pts;
    ac.initialize_container();
    ASSERT_EQ(ac.num_actuator_points(), num_pts);

    auto& pvec = data.position;
    // Inside
    pvec[0] = vs::Vector{64.0_rt, 64.0_rt, 64.0_rt};
    // On the low wall: inside
    pvec[1] = vs::Vector{64.0_rt, 64.0_rt, 0.0_rt};
    // Beyond a periodic boundary: inside
    pvec[2] = vs::Vector{130.0_rt, 64.0_rt, 64.0_rt};
    // Above the high wall
    pvec[3] = vs::Vector{64.0_rt, 64.0_rt, 130.0_rt};
    // On the high wall: outside, as for the particles
    pvec[4] = vs::Vector{64.0_rt, 64.0_rt, 128.0_rt};
    // Below the low wall
    pvec[5] = vs::Vector{64.0_rt, 64.0_rt, -1.0_rt};
    EXPECT_EQ(ac.num_points_outside_domain(), 3);

    // With the three points moved into the domain, every point gets the
    // velocity of the field
    pvec[3].z() = 100.0_rt;
    pvec[4].z() = 127.0_rt;
    pvec[5].z() = 1.0_rt;
    EXPECT_EQ(ac.num_points_outside_domain(), 0);
    ac.update_positions();
    ac.sample_fields(vel, density);
    const auto& vvec = data.velocity;
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e2_rt;
    for (int ip = 0; ip < num_pts; ++ip) {
        EXPECT_NEAR(vvec[ip].x(), 2.0_rt, tol) << "point " << ip;
    }
}

} // namespace kynema_sgf_tests
