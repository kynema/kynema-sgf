#include <numbers>
#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/turbulence/TurbulenceModel.H"
#include "src/fvm/stencils.H"
#include "src/utilities/tagging/CartBoxRefinement.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

// Velocity gradients: linear shear in z plus solid-body rotation about z
constexpr amrex::Real shear_x{0.5_rt};
constexpr amrex::Real shear_y{0.3_rt};
constexpr amrex::Real omega{0.2_rt};

/** Fill velocity, ghosts included, with a field of constant gradients
 *
 *  Ghosts beyond the non-periodic z faces hold the face value, which is what
 *  the one-sided wall stencils expect, so derivatives are exact everywhere.
 *
 *  \param vel [in,out] Velocity field
 */
void init_constant_gradient_velocity(kynema_sgf::Field& vel)
{
    const auto& mesh = vel.repo().mesh();
    for (int lev = 0; lev < vel.repo().num_active_levels(); ++lev) {
        const auto& geom = mesh.Geom(lev);
        const auto& dx = geom.CellSizeArray();
        const auto& problo = geom.ProbLoArray();
        const auto& probhi = geom.ProbHiArray();
        const auto& vel_arrs = vel(lev).arrays();
        amrex::ParallelFor(
            vel(lev), vel.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                const amrex::Real z = amrex::Clamp(
                    problo[2] + ((k + 0.5_rt) * dx[2]), problo[2], probhi[2]);
                vel_arrs[nbx](i, j, k, 0) = (shear_x * z) - (omega * y);
                vel_arrs[nbx](i, j, k, 1) = (shear_y * z) + (omega * x);
                vel_arrs[nbx](i, j, k, 2) = 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

/** Largest error in the Kosovic divNij on one level, over all ranks
 *
 *  Nij is uniform, so its divergence is zero away from the non-periodic z
 *  faces. In the cells next to those faces the one-sided stencils read a zero
 *  ghost, so the expected value is the difference between Nij = 0 on the face
 *  and the uniform Nij inside.
 *
 *  \param divNij [in] Scaled divergence of Nij computed by the model
 *  \param Nij [in] Non-linear stress tensor computed by the model
 *  \param lev [in] Level index
 *  \param scale [in] Model factor multiplying the raw divergence
 *  \return Maximum error in interior cells and in wall-adjacent cells
 */
amrex::Array<amrex::Real, 2> divnij_errors(
    const kynema_sgf::Field& divNij,
    const kynema_sgf::Field& Nij,
    const int lev,
    const amrex::Real scale)
{
    namespace stencil = kynema_sgf::fvm::stencil;
    const auto& geom = divNij.repo().mesh().Geom(lev);
    const amrex::Real idz = geom.InvCellSize(2);
    const int klo = geom.Domain().smallEnd(2);
    const int khi = geom.Domain().bigEnd(2);
    const amrex::Real lo_coeff =
        stencil::StencilKLO::c20 + stencil::StencilKLO::c21;
    const amrex::Real hi_coeff =
        stencil::StencilKHI::c21 + stencil::StencilKHI::c22;
    const auto& div_arrs = divNij(lev).const_arrays();
    const auto& nij_arrs = Nij(lev).const_arrays();
    auto errs = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMax, amrex::ReduceOpMax>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, divNij(lev),
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
            -> amrex::GpuTuple<amrex::Real, amrex::Real> {
            const bool at_wall = (k == klo) || (k == khi);
            const amrex::Real coeff =
                (k == klo) ? lo_coeff : ((k == khi) ? hi_coeff : 0.0_rt);
            amrex::Real err = 0.0_rt;
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                const amrex::Real expected =
                    scale * coeff *
                    nij_arrs[nbx](i, j, k, (n * AMREX_SPACEDIM) + 2) * idz;
                err = amrex::max<amrex::Real>(
                    err, std::abs(div_arrs[nbx](i, j, k, n) - expected));
            }
            if (at_wall) {
                return {0.0_rt, err};
            }
            return {err, 0.0_rt};
        });
    amrex::Array<amrex::Real, 2> out{amrex::get<0>(errs), amrex::get<1>(errs)};
    amrex::ParallelDescriptor::ReduceRealMax(out.data(), 2);
    return out;
}

/** Largest |Nij| in ghost cells beyond the non-periodic z faces, over all ranks
 *
 *  \param Nij [in] Non-linear stress tensor computed by the model
 *  \param lev [in] Level index
 *  \return Maximum magnitude of the wall ghost values
 */
amrex::Real max_wall_ghost_nij(const kynema_sgf::Field& Nij, const int lev)
{
    const auto& geom = Nij.repo().mesh().Geom(lev);
    const int klo = geom.Domain().smallEnd(2);
    const int khi = geom.Domain().bigEnd(2);
    const int ncomp = Nij.num_comp();
    const auto& nij_arrs = Nij(lev).const_arrays();
    amrex::Real val = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMax>{}, amrex::TypeList<amrex::Real>{},
        Nij(lev), amrex::IntVect(0, 0, 1),
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int k) -> amrex::GpuTuple<amrex::Real> {
            amrex::Real ghost = 0.0_rt;
            if ((k < klo) || (k > khi)) {
                for (int n = 0; n < ncomp; ++n) {
                    ghost = amrex::max<amrex::Real>(
                        ghost, std::abs(nij_arrs[nbx](i, j, k, n)));
                }
            }
            return ghost;
        });
    amrex::ParallelDescriptor::ReduceRealMax(val);
    return val;
}

} // namespace

/** Kosovic model on a mesh split into several boxes, periodic in x and y
 */
class KosovicNijTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_ncell, m_ncell, m_ncell}};
            pp.addarr("n_cell", ncell);
            pp.add("max_grid_size", m_max_grid_size);
            pp.add("blocking_factor", 4);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{16.0_rt, 16.0_rt, 16.0_rt}};
            amrex::Vector<int> periodic{{1, 1, 0}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("zlo");
            pp.add("type", std::string("slip_wall"));
        }
        {
            amrex::ParmParse pp("zhi");
            pp.add("type", std::string("slip_wall"));
        }
        {
            amrex::ParmParse pp("turbulence");
            pp.add("model", std::string("Kosovic"));
        }
        {
            amrex::ParmParse pp("Kosovic");
            pp.add("Cb", m_Cb);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"ABL"};
            pp.addarr("physics", physics);
            pp.add("density", m_rho);
            amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
            pp.addarr("velocity", vvec);
            amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -9.81_rt};
            pp.addarr("gravity", gvec);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("viscosity", 1.0e-5_rt);
            pp.add("reference_temperature", 300.0_rt);
        }
        {
            amrex::ParmParse pp("ABL");
            amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 4000.0_rt};
            pp.addarr("temperature_heights", t_hts);
            amrex::Vector<amrex::Real> t_vals{300.0_rt, 300.0_rt, 300.0_rt};
            pp.addarr("temperature_values", t_vals);
        }
    }

    //! Set up the model, impose the velocity and update the non-linear term
    void update_kosovic()
    {
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().init_physics();
        sim().create_turbulence_model();
        sim().turbulence_model().post_init_actions();

        init_constant_gradient_velocity(sim().repo().get_field("velocity"));
        sim().repo().get_field("density").setVal(m_rho);
        sim().repo().get_field("temperature").setVal(300.0_rt);

        sim().turbulence_model().update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    }

    //! Check divNij in the interior and next to the walls on every level
    void check_divnij()
    {
        const auto& repo = sim().repo();
        const auto& Nij = repo.get_field("Nij");
        const auto& divNij = repo.get_field("divNij");
        const amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

        // Nij is uniform, with N33 = 3 (du/dz^2 + dv/dz^2)
        for (int n = 0; n < Nij.num_comp(); ++n) {
            SCOPED_TRACE("Nij component " + std::to_string(n));
            EXPECT_NEAR(
                utils::field_min(Nij, n), utils::field_max(Nij, n), tol);
        }
        EXPECT_NEAR(
            utils::field_max(Nij, 8),
            3.0_rt * ((shear_x * shear_x) + (shear_y * shear_y)), tol);

        // Factor applied to the raw divergence in the Kosovic kernel
        const amrex::Real Cs_sqr = 8.0_rt * (1.0_rt + m_Cb) /
                                   (27.0_rt * std::numbers::pi_v<amrex::Real> *
                                    std::numbers::pi_v<amrex::Real>);
        const amrex::Real C1 =
            std::sqrt(960.0_rt) * m_Cb / (7.0_rt * (1.0_rt + m_Cb) * m_Sk);

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            SCOPED_TRACE("level " + std::to_string(lev));
            const auto& dx = mesh().Geom(lev).CellSizeArray();
            const amrex::Real ds = std::cbrt(dx[0] * dx[1] * dx[2]);
            const amrex::Real scale = m_rho * Cs_sqr * ds * ds * 0.25_rt * C1;

            // Interior cells, then cells next to the non-periodic z faces
            const auto errs = divnij_errors(divNij, Nij, lev, scale);
            EXPECT_LE(errs[0], tol);
            EXPECT_LE(errs[1], tol);
            // Nij stays zero beyond the non-periodic z faces
            EXPECT_EQ(max_wall_ghost_nij(Nij, lev), 0.0_rt);
        }
    }

    const int m_ncell{16};
    const int m_max_grid_size{8};
    const amrex::Real m_rho{1.2_rt};
    const amrex::Real m_Cb{0.36_rt};
    // Kosovic default backscatter parameter
    const amrex::Real m_Sk{0.5_rt};
};

/** Same setup with a refined patch touching the lower wall
 */
class KosovicNijAMRTest : public KosovicNijTest
{
protected:
    void populate_parameters() override
    {
        KosovicNijTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            pp.add("max_level", 1);
        }

        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "4.1 4.1 0.1 11.9 11.9 7.9" << '\n';

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<kynema_sgf::CartBoxRefinement> box_refine(
            new kynema_sgf::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }
};

TEST_F(KosovicNijTest, uniform_nij_divergence_across_boxes)
{
    update_kosovic();
    ASSERT_GT(mesh().boxArray(0).size(), 1);
    check_divnij();
}

TEST_F(KosovicNijAMRTest, uniform_nij_divergence_across_levels)
{
    update_kosovic();
    ASSERT_EQ(mesh().finestLevel(), 1);
    ASSERT_GT(mesh().boxArray(1).size(), 1);
    check_divnij();
}

} // namespace kynema_sgf_tests
