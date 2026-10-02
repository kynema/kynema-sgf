#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>

using namespace amrex::literals;

namespace {
// 100 m plateau for x in [449, 576], flat ground elsewhere
void write_terrain(const std::string& fname)
{
    std::ofstream os(fname);
    os << std::setprecision(17);
    os << "6\n2\n";
    os << "0.0\n448.0\n449.0\n576.0\n577.0\n1024.0\n";
    os << "0.0\n1024.0\n";
    os << "0.0\n0.0\n0.0\n0.0\n";
    os << "100.0\n100.0\n100.0\n100.0\n";
    os << "0.0\n0.0\n0.0\n0.0\n";
}
} // namespace

namespace kynema_sgf_tests {

// With the implicit immersed drag the factor 1 + C dt holds the velocity of the
// body at rest in the implicit diffusion solve. A scalar must not see it: the
// same factor would relax the scalar toward zero inside the terrain.
class ImmersedScalarDiffusionTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{1024.0_rt, 1024.0_rt, 512.0_rt}};
            pp.addarr("prob_hi", probhi);
            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("use_godunov", 1);
            pp.add("diffusion_type", 2);
            pp.add("density", m_rho);
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", m_dt);
        }
        {
            amrex::ParmParse pp("ImmersedTerrain");
            pp.add("implicit_projection", 1);
        }
    }
    const amrex::Real m_rho{1.2_rt};
    const amrex::Real m_dt{0.5_rt};
    const amrex::Real m_theta{300.0_rt};
};

TEST_F(ImmersedScalarDiffusionTest, scalar_in_body_is_not_damped)
{
    write_terrain("terrain.amrwind");
    populate_parameters();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().create_turbulence_model();
    sim().init_physics();
    auto& t_eqn = pde_mgr.register_transport_pde("Temperature");

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    auto terrain = std::make_unique<Terrain>(sim());
    terrain->initialize_fields(0, sim().repo().mesh().Geom(0));
    // The drag rate the projections and the velocity solve use is in place
    ASSERT_TRUE(sim().repo().field_exists("terrain_drag_rate"));
    const auto& rate = sim().repo().get_field("terrain_drag_rate");
    EXPECT_GT(utils::field_probe(rate, 0, 15, 10, 1), 0.0_rt);

    auto& density = sim().repo().get_field("density");
    density.setVal(m_rho);
    density.state(kynema_sgf::FieldState::Old).setVal(m_rho);
    density.state(kynema_sgf::FieldState::NPH).setVal(m_rho);
    auto& theta = t_eqn.fields().field;
    theta.setVal(m_theta);
    theta.state(kynema_sgf::FieldState::Old).setVal(m_theta);
    t_eqn.fields().mueff.setVal(0.1_rt);
    t_eqn.initialize();

    // A uniform scalar with no sources must come out of the implicit solve
    // unchanged, inside the terrain as well as in the fluid
    t_eqn.solve(m_dt);

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;
    EXPECT_NEAR(utils::field_probe(theta, 0, 15, 10, 1), m_theta, tol);
    EXPECT_NEAR(utils::field_probe(theta, 0, 15, 10, 3), m_theta, tol);
    EXPECT_NEAR(utils::field_probe(theta, 0, 5, 5, 8), m_theta, tol);
}

} // namespace kynema_sgf_tests
