#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/pp_utils.H"

#include "amr-wind/wind_energy/actuator/turbine/external/ExtTurbIface.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/fast_types.H"

#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

#define AW_ENABLE_OPENFAST_UTEST 0

namespace amr_wind_tests {
namespace {
class FastIfaceTest : public MeshTest
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
            amrex::Vector<amrex::Real> probhi{{128.0_rt, 128.0_rt, 256.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

} // namespace

TEST_F(FastIfaceTest, fast_init)
{
    initialize_mesh();
    pp_utils::default_time_inputs();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.0625_rt);
    }
    sim().time().parse_parameters();

    const int iproc = amrex::ParallelDescriptor::MyProc();
    ::ext_turb::FastTurbine fi;
    fi.tlabel = "T001";
    fi.tid_local = -1;
    fi.tid_global = iproc;
    fi.num_pts_blade = 5;
    fi.num_pts_tower = 5;
    fi.base_pos[0] = 64.0F;
    fi.base_pos[1] = 64.0F;
    fi.base_pos[2] = 0.0F;
    fi.input_file = "./fast_inp/nrel5mw.fst";
    fi.dt_cfd = 0.0625_rt;
    fi.start_time = 0.0_rt;
    fi.stop_time = 0.625_rt;

    ::ext_turb::ExtTurbIface<
        ::ext_turb::FastTurbine, ::ext_turb::FastSolverData>
        fast(sim());
    fast.parse_inputs(sim(), "OpenFAST");
    fast.register_turbine(fi);

    EXPECT_EQ(fi.tid_local, 0);

#if AW_ENABLE_OPENFAST_UTEST
    fast.init_turbine(fi.tid_local);
    EXPECT_NEAR(
        fi.dt_ext, 0.00625_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    EXPECT_EQ(fi.num_substeps, 10);
    EXPECT_EQ(fi.num_blades, 3);
    EXPECT_TRUE(fi.is_solution0);

    std::fill(fi.from_cfd.u, fi.from_cfd.u + fi.from_cfd.u_Len, 6.0_rt);
    std::fill(fi.from_cfd.v, fi.from_cfd.v + fi.from_cfd.v_Len, 0.0_rt);
    std::fill(fi.from_cfd.w, fi.from_cfd.w + fi.from_cfd.w_Len, 0.0_rt);
    fast.init_solution(fi.tid_local);
    EXPECT_FALSE(fi.is_solution0);

    for (int i = 0; i < 2; ++i) {
        fast.advance_turbine(fi.tid_local);
    }
    EXPECT_EQ(fi.time_index, 20);
#else
    GTEST_SKIP();
#endif
}

TEST_F(FastIfaceTest, fast_replay)
{
    initialize_mesh();
    pp_utils::default_time_inputs();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.0625_rt);
    }
    sim().time().parse_parameters();

    const int iproc = amrex::ParallelDescriptor::MyProc();
    ::ext_turb::FastTurbine fi;
    fi.tlabel = "T001";
    fi.tid_local = -1;
    fi.tid_global = iproc;
    fi.num_pts_blade = 5;
    fi.num_pts_tower = 5;
    fi.base_pos[0] = 64.0F;
    fi.base_pos[1] = 64.0F;
    fi.base_pos[2] = 0.0F;
    fi.input_file = "./fast_inp/nrel5mw1.fst";
    fi.start_time = 0.125_rt;
    fi.stop_time = 0.625_rt;
    fi.dt_cfd = 0.0625_rt;
    fi.sim_mode = ::ext_turb::SimMode::replay;

    ::ext_turb::ExtTurbIface<
        ::ext_turb::FastTurbine, ::ext_turb::FastSolverData>
        fast(sim());
    fast.parse_inputs(sim(), "OpenFAST");
    fast.register_turbine(fi);

    EXPECT_EQ(fi.tid_local, 0);
    EXPECT_TRUE(fi.is_solution0);

#if AW_ENABLE_OPENFAST_UTEST
    fast.init_turbine(fi.tid_local);
    EXPECT_EQ(fi.time_index, 20);
    EXPECT_FALSE(fi.is_solution0);

    fast.advance_turbine(fi.tid_local);
    EXPECT_EQ(fi.time_index, 30);
#else
    GTEST_SKIP();
#endif
}

} // namespace amr_wind_tests
