#include "aw_test_utils/MeshTest.H"
#include "amr-wind/equation_systems/icns/icns_advection.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class ICNSConstDensTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        // Specify uniform density according to how it is done in input files
        {
            amrex::ParmParse pp("incflo");
            pp.add("density", m_rho_0);
        }
    }

    void testing_density()
    {
        initialize_mesh();

        // Initialize MAC projection operator
        const auto& mco = amr_wind::pde::MacProjOp(
            sim().repo(), sim().physics_manager(), false, false, false, false);
        // Get background density and check
        const amrex::Real rho0 = mco.rho0();
        EXPECT_EQ(rho0, m_rho_0);
    }
    const amrex::Real m_rho_0 = 2.0_rt;
};

TEST_F(ICNSConstDensTest, nonunity) { testing_density(); }

} // namespace amr_wind_tests
