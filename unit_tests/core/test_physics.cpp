#include "ks_test_utils/MeshTest.H"
#include "src/core/Physics.H"
#include "physics_test_utils.H"

namespace kynema_sgf_tests {

class PhysicsTest : public MeshTest
{};

TEST_F(PhysicsTest, physics_example)
{
    PhysicsEx obj(sim());
    EXPECT_EQ("PhysicsEx", obj.identifier());
}

TEST_F(PhysicsTest, physics_interface)
{
    initialize_mesh();
    auto& phy_mgr = sim().physics_manager();

    EXPECT_FALSE(phy_mgr.contains(PhysicsEx::identifier()));
    phy_mgr.create("PhysicsEx", sim());
    EXPECT_TRUE(phy_mgr.contains(PhysicsEx::identifier()));

    if (phy_mgr.contains(PhysicsEx::identifier())) {
        const auto& pex = phy_mgr.get<PhysicsEx>();
        EXPECT_NE(&pex, nullptr);
    }
}

TEST_F(PhysicsTest, physics_init_inputs)
{
    initialize_mesh();
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"PhysicsEx"};
        pp.addarr("physics", physics);
    }

    auto& phy_mgr = sim().physics_manager();
    EXPECT_FALSE(phy_mgr.contains(PhysicsEx::identifier()));
    sim().init_physics();
    EXPECT_TRUE(phy_mgr.contains(PhysicsEx::identifier()));

    if (phy_mgr.contains(PhysicsEx::identifier())) {
        const auto& pex = phy_mgr.get<PhysicsEx>();
        EXPECT_NE(&pex, nullptr);
    }
}

TEST_F(PhysicsTest, physics_test_duplicates)
{
    initialize_mesh();
    auto& phy_mgr = sim().physics_manager();

    EXPECT_FALSE(phy_mgr.contains(PhysicsEx::identifier()));
    auto& pex1 = phy_mgr.create("PhysicsEx", sim());
    auto& pex2 = phy_mgr.create("PhysicsEx", sim());
    EXPECT_EQ(phy_mgr.objects().size(), 1U);
    EXPECT_TRUE(phy_mgr.contains(PhysicsEx::identifier()));

    EXPECT_EQ(&pex1, &pex2);
}

} // namespace kynema_sgf_tests
