#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/incflo.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

using namespace amrex::literals;

namespace {
//! Flat terrain of the given height over the unit square
void write_flat_terrain(const std::string& fname, const double height)
{
    std::ofstream os(fname);
    os << std::setprecision(17);
    os << "2\n2\n";
    os << "0.0\n1.0\n";
    os << "0.0\n1.0\n";
    for (int n = 0; n < 4; ++n) {
        os << height << "\n";
    }
}

/** Nodal projection of a uniform horizontal velocity with the implicit
 *  immersed drag. The field depends on z only, so it is divergence free and
 *  the projection leaves it unchanged except for the implicit drag factor.
 *
 *  \param dt Time step of the projection
 *  \param dt_nm1 Previous time step; below 10 dt_nm1 the projection uses its
 *  small-dt (delta) form
 *  \param u_body Returned velocity in a cell inside the terrain
 *  \param u_fluid Returned velocity in a fluid cell
 */
void project_uniform_flow(
    const amrex::Real dt,
    const amrex::Real dt_nm1,
    amrex::Real& u_body,
    amrex::Real& u_fluid)
{
    incflo my_incflo;
    my_incflo.init_mesh();
    auto& repo = my_incflo.sim().repo();
    auto& density = repo.get_field("density");
    auto& velocity = repo.get_field("velocity");
    auto& velocity_old = velocity.state(kynema_sgf::FieldState::Old);
    density.setVal(1.0_rt);
    repo.get_field("gp").setVal(0.0_rt);
    repo.get_field("p").setVal(0.0_rt);
    // u** = u^n = (1, 0, 0) everywhere
    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        velocity(lev).setVal(0.0_rt);
        velocity(lev).setVal(1.0_rt, 0, 1, velocity.num_grow());
        velocity_old(lev).setVal(0.0_rt);
        velocity_old(lev).setVal(1.0_rt, 0, 1, velocity_old.num_grow());
    }

    auto& time = my_incflo.sim().time();
    time.delta_t() = dt;
    time.delta_t_nm1() = dt_nm1;
    my_incflo.ApplyProjection(density.vec_const_ptrs(), 1.0_rt, dt, false);

    u_body = kynema_sgf_tests::utils::field_probe(velocity, 0, 1, 1, 2, 0);
    u_fluid = kynema_sgf_tests::utils::field_probe(velocity, 0, 1, 1, 12, 0);
}
} // namespace

namespace kynema_sgf_tests {

// The implicit immersed drag divides the predicted velocity by 1 + C dt in the
// nodal projection, consistently with the effective density rho (1 + C dt) of
// the Poisson coefficient. This holds in the full form and in the small-dt
// form that the projection switches to when the time step drops sharply.
class ImmersedProjectionTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{1.0_rt, 1.0_rt, 1.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            amrex::Vector<int> periodic{{1, 1, 0}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("use_godunov", 1);
            amrex::Vector<std::string> physics{"ImmersedTerrain"};
            pp.addarr("physics", physics);
        }
        {
            amrex::ParmParse pp("ImmersedTerrain");
            pp.add("terrain_file", m_terrain_file);
            pp.add("drag_weight", std::string("center"));
            pp.add("implicit_projection", 1);
        }
        {
            amrex::ParmParse pp("ImmersedDragForcing");
            pp.add("drag_coefficient", m_cd);
        }
        amrex::ParmParse ppzlo("zlo");
        ppzlo.add("type", std::string("slip_wall"));
        amrex::ParmParse ppzhi("zhi");
        ppzhi.add("type", std::string("pressure_outflow"));
    }

    // dz = 1/16; the terrain at 0.3 fills k = 0..3 and 80 % of k = 4, all
    // solid with center weighting
    const std::string m_terrain_file{"terrain_projection.amrwind"};
    const amrex::Real m_height{0.3_rt};
    const amrex::Real m_cd{10.0_rt};
    const amrex::Real m_dt{1.0e-3_rt};
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

TEST_F(ImmersedProjectionTest, full_form)
{
    write_flat_terrain(m_terrain_file, m_height);
    populate_parameters();
    initialize_mesh();
    amrex::Real u_body = 0.0_rt;
    amrex::Real u_fluid = 0.0_rt;
    project_uniform_flow(m_dt, m_dt, u_body, u_fluid);
    const amrex::Real rate = m_cd * 16.0_rt;
    EXPECT_NEAR(u_body, 1.0_rt / (1.0_rt + (rate * m_dt)), m_tol);
    EXPECT_NEAR(u_fluid, 1.0_rt, m_tol);
}

TEST_F(ImmersedProjectionTest, small_dt_form)
{
    write_flat_terrain(m_terrain_file, m_height);
    populate_parameters();
    initialize_mesh();
    amrex::Real u_body = 0.0_rt;
    amrex::Real u_fluid = 0.0_rt;
    // dt below a tenth of the previous step: projection of u** - u^n
    project_uniform_flow(m_dt, 100.0_rt * m_dt, u_body, u_fluid);
    const amrex::Real rate = m_cd * 16.0_rt;
    EXPECT_NEAR(u_body, 1.0_rt / (1.0_rt + (rate * m_dt)), m_tol);
    EXPECT_NEAR(u_fluid, 1.0_rt, m_tol);
}

} // namespace kynema_sgf_tests
