#include <array>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "ks_test_utils/AmrexTest.H"
#include "test_ext_turbine_utils.H"
#include "src/actuator/turbine/external/turbine_external_ops.H"
#include "src/actuator/turbine/ActSrcDiskOp_Turbine.H"
#include "AMReX_Exception.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

// Stand-ins for the turbine and solver data of an external solver. The
// turbine holds the buffers exchanged with the solver: the fluid velocity
// at the velocity points and the velocity of the force points.
struct MockExtTurbine
{
    float* fluid_velocity(const int dir) { return fluid_vel[dir].data(); }
    float* solid_velocity(const int dir) { return solid_vel[dir].data(); }
    [[nodiscard]] int length_fluid_velocity(const int dir) const
    {
        return static_cast<int>(fluid_vel[dir].size());
    }

    std::array<std::vector<float>, 3> fluid_vel;
    std::array<std::vector<float>, 3> solid_vel;
};

struct MockExtSolverData
{};

using MockViewData = ExtTurbineViewData<MockExtTurbine, MockExtSolverData>;

// A solver whose blade points have the local frame of Kynema-FMB
struct MockSpanFirstTurbine : public MockExtTurbine
{
    static constexpr ::ext_turb::BladeFrame blade_frame{
        ::ext_turb::BladeFrame::span_chord_normal};
};

using MockSpanFirstViewData =
    ExtTurbineViewData<MockSpanFirstTurbine, MockExtSolverData>;

//! Inputs (chord, thickness, span) = (1, 2, 3): check every blade and
//! tower point in the given local frame
template <typename ViewData>
void check_epsilon_frame(
    ViewData& data,
    const amrex::Real along_x,
    const amrex::Real along_y,
    const amrex::Real along_z)
{
    namespace vs = kynema_sgf::vs;
    init_hub_first_points(data, 3, 5, 4);
    kynema_sgf::actuator::external::make_component_views(data);
    auto& meta = data.meta();
    meta.eps_inp = vs::Vector{1.0_rt, 2.0_rt, 3.0_rt};
    meta.eps_tower = vs::Vector{1.0_rt, 2.0_rt, 3.0_rt};
    kynema_sgf::actuator::external::init_epsilon(data);

    const vs::Vector expected{along_x, along_y, along_z};
    for (int ib = 0; ib < meta.num_blades; ++ib) {
        for (const auto& eps : meta.blades[ib].epsilon) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                EXPECT_EQ(eps[n], expected[n]) << "blade " << ib;
            }
        }
    }
    for (const auto& eps : meta.tower.epsilon) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_EQ(eps[n], expected[n]) << "tower";
        }
    }
}

} // namespace

TEST(TurbineExternalViews, hub_first_velocity_points)
{
    MockViewData data;
    init_hub_first_points(data, 3, 5, 5);

    // The hub is first in the velocity points: the views skip it, as in the
    // force points
    kynema_sgf::actuator::external::make_component_views(data);

    check_blade_views_aligned(data);
}

TEST(TurbineExternalViews, fllc_relative_velocity_with_different_point_counts)
{
    namespace vs = kynema_sgf::vs;
    const int num_blades = 3;
    const int npb = 4; // force points on a blade
    const int nvb = 6; // velocity points on a blade
    MockViewData data;
    init_hub_first_points(data, num_blades, npb, nvb, 2);
    auto& grid = data.grid();
    auto& meta = data.meta();
    kynema_sgf::actuator::external::make_component_views(data);
    data.info().is_root_proc = true;

    // Fluid velocity last sent to the external solver: linear in space, so
    // that its interpolation along a blade is exact. Velocity of the force
    // points: different for every point. Velocity sampled now: not linear.
    const auto fluid_u = [](const vs::Vector& x) {
        return 1.0_rt + (2.0_rt * x.y()) + (3.0_rt * x.z());
    };
    const int npts = static_cast<int>(grid.pos.size());
    const int npts_vel = static_cast<int>(grid.vel_pos.size());
    auto& ext = meta.ext_data;
    for (int dir = 0; dir < 3; ++dir) {
        ext.fluid_vel[dir].assign(npts_vel, 0.0F);
        ext.solid_vel[dir].assign(npts, 0.0F);
    }
    for (int i = 0; i < npts_vel; ++i) {
        ext.fluid_vel[0][i] = static_cast<float>(fluid_u(grid.vel_pos[i]));
        grid.vel[i] = vs::Vector{
            static_cast<amrex::Real>(7 + ((i * i) % 5)), 0.0_rt, 0.0_rt};
    }
    for (int i = 0; i < npts; ++i) {
        ext.solid_vel[0][i] = 0.25F * static_cast<float>(i);
    }

    // Correction initialized and no force on the blades
    meta.fllc.resize(num_blades);
    for (int ib = 0; ib < num_blades; ++ib) {
        meta.fllc[ib].epsilon = 2.0_rt;
        kynema_sgf::actuator::fllc_init(meta.fllc[ib], meta.blades[ib], 1.0_rt);
    }

    kynema_sgf::actuator::ops::update_fields_op(data);

    // Relative velocity at force point k of every blade: the fluid velocity
    // at that point minus the velocity of that point
    const amrex::Real tol = 1.0e-5_rt;
    for (int ib = 0; ib < num_blades; ++ib) {
        for (int k = 0; k < npb; ++k) {
            const int ip = 1 + (ib * npb) + k;
            const auto gold = fluid_u(grid.pos[ip]) -
                              (0.25_rt * static_cast<amrex::Real>(ip));
            EXPECT_NEAR(meta.blades[ib].vel_rel[k].x(), gold, tol)
                << "blade " << ib << ", force node " << k;
            EXPECT_EQ(meta.blades[ib].vel_rel[k].y(), 0.0_rt);
            EXPECT_EQ(meta.blades[ib].vel_rel[k].z(), 0.0_rt);
        }
    }
    // Without force there is no correction: the velocity sent to the
    // external solver is the sampled one, at every velocity point
    for (int i = 0; i < npts_vel; ++i) {
        EXPECT_EQ(ext.fluid_vel[0][i], static_cast<float>(7 + ((i * i) % 5)))
            << "velocity point " << i;
        EXPECT_EQ(ext.fluid_vel[1][i], 0.0F);
        EXPECT_EQ(ext.fluid_vel[2][i], 0.0F);
    }
}

TEST(TurbineExternalViews, epsilon_in_frame_of_openfast)
{
    // Local axes: normal to the chord, along the chord, along the span
    MockViewData data;
    check_epsilon_frame(data, 2.0_rt, 1.0_rt, 3.0_rt);
}

TEST(TurbineExternalViews, epsilon_in_frame_of_kynema_fmb)
{
    // Local axes: along the span, along the chord, normal to the chord
    MockSpanFirstViewData data;
    check_epsilon_frame(data, 3.0_rt, 1.0_rt, 2.0_rt);
}

TEST(TurbineExternalViews, epsilon_component_of_each_input)
{
    using ::ext_turb::BladeFrame;
    using kynema_sgf::actuator::external::epsilon_component;
    // Inputs: chord, thickness, span
    EXPECT_EQ(epsilon_component(BladeFrame::normal_chord_span, 0), 1);
    EXPECT_EQ(epsilon_component(BladeFrame::normal_chord_span, 1), 0);
    EXPECT_EQ(epsilon_component(BladeFrame::normal_chord_span, 2), 2);
    EXPECT_EQ(epsilon_component(BladeFrame::span_chord_normal, 0), 1);
    EXPECT_EQ(epsilon_component(BladeFrame::span_chord_normal, 1), 2);
    EXPECT_EQ(epsilon_component(BladeFrame::span_chord_normal, 2), 0);
    static_assert(
        ::ext_turb::blade_frame_of<MockExtTurbine>() ==
        BladeFrame::normal_chord_span);
    static_assert(
        ::ext_turb::blade_frame_of<MockSpanFirstTurbine>() ==
        BladeFrame::span_chord_normal);
}

TEST(TurbineExternalViews, epsilon_must_be_positive_on_blades)
{
    MockViewData data;
    init_hub_first_points(data, 3, 5, 0);
    kynema_sgf::actuator::external::make_component_views(data);

    // None of the epsilon inputs given: the Gaussian of every blade point
    // has no width, and its force would not reach the mesh
    EXPECT_THROW(
        kynema_sgf::actuator::external::init_epsilon(data),
        amrex::RuntimeError);
}

TEST(TurbineExternalViews, epsilon_must_be_positive_on_tower)
{
    namespace vs = kynema_sgf::vs;
    MockViewData data;
    init_hub_first_points(data, 3, 5, 4);
    kynema_sgf::actuator::external::make_component_views(data);

    // Tower points, and an epsilon for the blades only
    data.meta().eps_chord = vs::Vector{0.25_rt, 0.25_rt, 0.25_rt};
    EXPECT_THROW(
        kynema_sgf::actuator::external::init_epsilon(data),
        amrex::RuntimeError);
}

TEST(TurbineExternalViews, epsilon_from_chord_and_tower_inputs)
{
    namespace vs = kynema_sgf::vs;
    const int num_blades = 3;
    const int npb = 5;
    const int npt = 4;
    MockViewData data;
    init_hub_first_points(data, num_blades, npb, npt);
    auto& meta = data.meta();
    kynema_sgf::actuator::external::make_component_views(data);

    meta.eps_chord = vs::Vector{0.25_rt, 0.25_rt, 0.25_rt};
    meta.eps_tower = vs::Vector{2.0_rt, 2.0_rt, 2.0_rt};
    kynema_sgf::actuator::external::init_epsilon(data);

    // Unit chord on the blades
    for (int ib = 0; ib < num_blades; ++ib) {
        for (int k = 0; k < npb; ++k) {
            EXPECT_EQ(meta.blades[ib].epsilon[k].x(), 0.25_rt);
        }
    }
    for (int j = 0; j < npt; ++j) {
        EXPECT_EQ(meta.tower.epsilon[j].x(), 2.0_rt);
    }
}

TEST(TurbineExternalViews, bound_box_holds_the_gaussians)
{
    namespace vs = kynema_sgf::vs;
    MockViewData data;
    init_hub_first_points(data, 3, 5, 4);
    auto& meta = data.meta();
    kynema_sgf::actuator::external::make_component_views(data);

    // Box: 1.25 radii around a rotor of radius 5 at the origin; points up to 5
    // from the hub
    data.info().bound_box =
        amrex::RealBox(-6.25_rt, -6.25_rt, -6.25_rt, 6.25_rt, 6.25_rt, 6.25_rt);
    meta.eps_inp = vs::Vector{2.0_rt, 2.0_rt, 2.0_rt};
    kynema_sgf::actuator::external::init_epsilon(data);

    // Gaussian cutoff at 4 epsilon: the box has to reach that far
    const amrex::Real reach = 4.0_rt * 2.0_rt;
    const auto& bbox = data.info().bound_box;
    for (const auto& pos : data.grid().pos) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_LE(bbox.lo(n), pos[n] - reach);
            EXPECT_GE(bbox.hi(n), pos[n] + reach);
        }
    }
}

//! Disk width normal to the rotor for inputs (chord, thickness, span) =
//! (1, 2, 3), in the solver frame
template <typename ViewData>
amrex::Real disk_normal_epsilon()
{
    using MetaType =
        std::remove_cvref_t<decltype(std::declval<ViewData&>().meta())>;
    using SolverTurbine =
        std::remove_cvref_t<decltype(std::declval<MetaType&>().ext_data)>;
    kynema_sgf::vs::Vector eps{1.0_rt, 2.0_rt, 3.0_rt};
    kynema_sgf::actuator::external::epsilon_to_frame(
        eps, ::ext_turb::blade_frame_of<SolverTurbine>());
    return eps
        [kynema_sgf::actuator::ops::disk_normal_epsilon_component<MetaType>()];
}

TEST(TurbineExternalViews, disk_width_is_thickness_value)
{
    // Thickness input for both frames
    EXPECT_EQ(disk_normal_epsilon<MockViewData>(), 2.0_rt);
    EXPECT_EQ(disk_normal_epsilon<MockSpanFirstViewData>(), 2.0_rt);

    // Turbines without an external solver: first component
    struct NoExtMeta
    {};
    static_assert(
        kynema_sgf::actuator::ops::disk_normal_epsilon_component<NoExtMeta>() ==
        0);
}

} // namespace kynema_sgf_tests
