#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"

namespace amr_wind_tests {

class PDETest : public MeshTest
{};

TEST_F(PDETest, test_pde_create_godunov)
{
    amrex::ParmParse pp("incflo");
    pp.add("use_godunov", 1);

    initialize_mesh();

    auto& pde_mgr = mesh().sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");

    EXPECT_EQ(pde_mgr.scalar_eqns().size(), 1);

    EXPECT_EQ(mesh().field_repo().num_fields(), 22);
}

TEST_F(PDETest, test_pde_create_mol)
{
    amrex::ParmParse pp("incflo");
    pp.add("use_godunov", 0);

    initialize_mesh();

    auto& pde_mgr = mesh().sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");

    EXPECT_EQ(pde_mgr.scalar_eqns().size(), 1);

    EXPECT_EQ(mesh().field_repo().num_fields(), 26);
}

} // namespace amr_wind_tests
