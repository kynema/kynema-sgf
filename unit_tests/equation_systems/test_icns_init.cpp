#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/utilities/tagging/CartBoxRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class ICNSInitTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{1.0_rt, 1.0_rt, 1.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_nx + 1, m_nx + 1, m_nx + 1}};
            pp.add("max_level", 1);
            pp.add("max_grid_size", m_nx + 1);
            pp.add("blocking_factor", 2);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        // Create the "input file"
        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "0.8 0.5 0.5 0.9 0.5 0.5" << '\n';

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
            new amr_wind::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    const int m_nx{3};
};

TEST_F(ICNSInitTest, 2level)
{
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.icns().initialize();
}

TEST_F(ICNSInitTest, generic_2level)
{
    populate_parameters();
    initialize_mesh();
    auto& repo = sim().repo();
    auto& generic_field = repo.declare_field("generic", 1, 1, 1);
    generic_field.set_default_fillpatch_bc(sim().time());
    generic_field.fillpatch(0.0_rt);
}

} // namespace amr_wind_tests
