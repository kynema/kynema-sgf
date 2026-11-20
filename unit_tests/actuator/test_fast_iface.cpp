#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/pp_utils.H"

#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "src/actuator/turbine/fast/fast_types.H"
#include "src/utilities/ncutils/nc_interface.H"
#include "test_ext_turbine_utils.H"

#include <algorithm>
#include <string>
#include <vector>
#include "AMReX_FileSystem.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

#define KS_ENABLE_OPENFAST_UTEST 0

namespace kynema_sgf_tests {
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

#ifdef KYNEMA_SGF_USE_NETCDF
//! Velocity file functions of the OpenFAST interface, open to the tests
class FastVelocityFile
    : public ::ext_turb::
          ExtTurbIface<::ext_turb::FastTurbine, ::ext_turb::FastSolverData>
{
public:
    using FastIface = ::ext_turb::
        ExtTurbIface<::ext_turb::FastTurbine, ::ext_turb::FastSolverData>;

    explicit FastVelocityFile(const kynema_sgf::CFDSim& sim) : FastIface(sim) {}

    using FastIface::prepare_netcdf_file;
    using FastIface::read_velocity_data;
    using FastIface::write_velocity_data;

    std::string file_name(const ::ext_turb::FastTurbine& fi) const
    {
        return m_solver_data.output_dir + "/" + fi.tlabel + ".nc";
    }
};

/** Turbine with its velocity arrays held here, not by OpenFAST
 *
 *  Two OpenFAST steps of 0.25 s per CFD step.
 */
struct VelocityFileTurbine
{
    explicit VelocityFileTurbine(const std::string& label)
    {
        fi.tlabel = label;
        fi.dt_ext = 0.25_rt;
        fi.num_substeps = 2;
        fi.to_cfd.pxVel = pos.data();
        fi.to_cfd.pyVel = pos.data();
        fi.to_cfd.pzVel = pos.data();
        fi.from_cfd.u = u.data();
        fi.from_cfd.v = v.data();
        fi.from_cfd.w = w.data();
        fi.from_cfd.u_Len = npts;
        fi.from_cfd.v_Len = npts;
        fi.from_cfd.w_Len = npts;
    }

    //! Record of a CFD step, one velocity at all points
    void write_step(FastVelocityFile& fast, const int step, const float vel)
    {
        fi.time_index = step * fi.num_substeps;
        std::fill(u.begin(), u.end(), vel);
        fast.write_velocity_data(fi);
    }

    static constexpr int npts = 4;
    std::vector<float> pos = std::vector<float>(npts, 0.0F);
    std::vector<float> u = std::vector<float>(npts, 0.0F);
    std::vector<float> v = std::vector<float>(npts, 0.0F);
    std::vector<float> w = std::vector<float>(npts, 0.0F);
    ::ext_turb::FastTurbine fi;
};
#endif

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

#if KS_ENABLE_OPENFAST_UTEST
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

#if KS_ENABLE_OPENFAST_UTEST
    fast.init_turbine(fi.tid_local);
    EXPECT_EQ(fi.time_index, 20);
    EXPECT_FALSE(fi.is_solution0);

    fast.advance_turbine(fi.tid_local);
    EXPECT_EQ(fi.time_index, 30);
#else
    GTEST_SKIP();
#endif
}

TEST_F(FastIfaceTest, fast_velocity_points)
{
    // OpenFAST lists the hub first in its force points and in its velocity
    // points. The blade views must pair every force point with its own
    // velocity point.
    ExtTurbineViewData<::ext_turb::FastTurbine, ::ext_turb::FastSolverData>
        data;
    init_hub_first_points(data, 3, 5, 5);
    ::kynema_sgf::actuator::external::make_component_views(data);

    check_blade_views_aligned(data);
}

TEST_F(FastIfaceTest, velocity_record_is_time_step_of_turbine)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    initialize_mesh();
    // Velocity file: named after the turbine label, written by the turbine
    // root only. The I/O rank plays the root.
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    FastVelocityFile fast(sim());
    VelocityFileTurbine turbine("utest_velocity_record");
    const std::string fname = fast.file_name(turbine.fi);
    // prepare_netcdf_file keeps an existing file
    amrex::FileSystem::Remove(fname);
    fast.prepare_netcdf_file(turbine.fi);

    // 6 CFD steps; velocity n at step n
    const int nsteps = 6;
    for (int step = 0; step < nsteps; ++step) {
        turbine.write_step(fast, step, static_cast<float>(step));
    }
    // Restart from step 3: steps 3 and 4 again, other velocities
    turbine.write_step(fast, 3, 30.0F);
    turbine.write_step(fast, 4, 40.0F);

    {
        const auto ncf = ncutils::NCFile::open(fname);
        // No record added for the repeated steps
        ASSERT_EQ(static_cast<int>(ncf.dim("num_time_steps").len()), nsteps);

        const int npts = VelocityFileTurbine::npts;
        std::vector<float> time(nsteps, -1.0F);
        std::vector<float> uvel(nsteps * npts, -1.0F);
        ncf.var("time").get(time.data(), {0}, {nsteps});
        ncf.var("uvel").get(uvel.data(), {0, 0}, {nsteps, npts});
        const std::vector<float> uvel_expected{0.0F,  1.0F,  2.0F,
                                               30.0F, 40.0F, 5.0F};
        for (int step = 0; step < nsteps; ++step) {
            EXPECT_FLOAT_EQ(time[step], 0.5F * static_cast<float>(step));
            EXPECT_FLOAT_EQ(uvel[step * npts], uvel_expected[step]);
        }
    }
    amrex::FileSystem::Remove(fname);
#else
    GTEST_SKIP();
#endif
}

TEST_F(FastIfaceTest, replay_rejects_velocity_record_of_another_time)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    initialize_mesh();
    // Velocity file: named after the turbine label, written by the turbine
    // root only. The I/O rank plays the root.
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    FastVelocityFile fast(sim());
    VelocityFileTurbine turbine("utest_velocity_replay");
    const std::string fname = fast.file_name(turbine.fi);
    amrex::FileSystem::Remove(fname);
    fast.prepare_netcdf_file(turbine.fi);

    // File of a run restarted from step 5 in a new directory: steps 5 to 7
    turbine.write_step(fast, 5, 5.0F);
    turbine.write_step(fast, 6, 6.0F);
    turbine.write_step(fast, 7, 7.0F);

    {
        const auto ncf = ncutils::NCFile::open(fname);
        // Step 6 is in the file
        FastVelocityFile::read_velocity_data(turbine.fi, ncf, 6);
        EXPECT_FLOAT_EQ(turbine.u[0], 6.0F);
        // Step 2 is not: abort
        EXPECT_THROW(
            FastVelocityFile::read_velocity_data(turbine.fi, ncf, 2),
            amrex::RuntimeError);
        // Other time step: step 6 at another time than the record
        turbine.fi.dt_ext = 0.5_rt;
        EXPECT_THROW(
            FastVelocityFile::read_velocity_data(turbine.fi, ncf, 6),
            amrex::RuntimeError);
    }
    amrex::FileSystem::Remove(fname);
#else
    GTEST_SKIP();
#endif
}

} // namespace kynema_sgf_tests
