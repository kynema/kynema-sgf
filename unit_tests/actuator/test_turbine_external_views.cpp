#include "ks_test_utils/AmrexTest.H"
#include "test_ext_turbine_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

// Stand-ins for the turbine and solver data of an external solver
struct MockExtTurbine
{};

struct MockExtSolverData
{};

using MockViewData = ExtTurbineViewData<MockExtTurbine, MockExtSolverData>;

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

} // namespace kynema_sgf_tests
