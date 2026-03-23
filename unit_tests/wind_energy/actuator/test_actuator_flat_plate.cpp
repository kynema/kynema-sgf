#include <algorithm>

#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/ActSrcLineOp.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/core/gpu_utils.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
namespace {
class ActFlatPlateTest : public MeshTest
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
            amrex::Vector<amrex::Real> probhi{{32.0_rt, 32.0_rt, 32.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace act = amr_wind::actuator;
namespace vs = amr_wind::vs;

struct FlatPlate : public act::WingType
{
    using InfoType = act::ActInfo;
    using GridType = act::ActGrid;
    using MetaType = act::WingBaseData;
    using DataType = act::ActDataHolder<FlatPlate>;

    static std::string identifier() { return "TestFlatPlate"; }
};

} // namespace

} // namespace amr_wind_tests

namespace amr_wind::actuator {

namespace ops {

template <>
struct UseDefaultOp<
    ::amr_wind_tests::FlatPlate,
    ::amr_wind::actuator::ActSrcLine>
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = false;
    static constexpr bool compute_force = false;
    static constexpr bool process_outputs = true;
};

template <typename SrcTrait>
struct ReadInputsOp<::amr_wind_tests::FlatPlate, SrcTrait>
{
    void operator()(
        ::amr_wind_tests::FlatPlate::DataType& data,
        const amr_wind::actuator::utils::ActParser& pp)
    {
        // Copied from flat_plate_ops, but with unit chord assumed
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.get_either("epsilon", wdata.eps_inp);
        pp.get("pitch", wdata.pitch);

        amrex::Real max_eps = *std::ranges::max_element(wdata.eps_inp);
        amrex::Real search_radius = max_eps * 3.0_rt;
        const auto& p1 = wdata.start;
        const auto& p2 = wdata.end;
        // clang-format off
        info.bound_box = amrex::RealBox(
            amrex::min<amrex::Real>(p1.x(), p2.x()) - search_radius,
            amrex::min<amrex::Real>(p1.y(), p2.y()) - search_radius,
            amrex::min<amrex::Real>(p1.z(), p2.z()) - search_radius,
            amrex::max<amrex::Real>(p1.x(), p2.x()) + search_radius,
            amrex::max<amrex::Real>(p1.y(), p2.y()) + search_radius,
            amrex::max(p1.z(), p2.z()) + search_radius
        );
        // clang-format on
    }
};

template <typename SrcTrait>
struct InitDataOp<::amr_wind_tests::FlatPlate, SrcTrait>
{
    void operator()(::amr_wind_tests::FlatPlate::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();

        amr_wind::actuator::wing::init_data_structures(wdata, grid);

        {
            // Check that the wing coordinates were initialized correctly
            auto wing_len = vs::mag(grid.pos.back() - grid.pos.front());
            auto angle1 =
                std::acos(grid.orientation[0].x().unit() & vs::Vector::ihat());
            auto angle2 = std::asin(
                (grid.orientation[0].x().unit() ^ vs::Vector::ihat()) &
                vs::Vector::jhat());
            auto angle_gold = ::amr_wind::utils::radians(wdata.pitch);
            EXPECT_NEAR(
                wing_len, 8.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                angle1, angle_gold,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                angle2, -angle_gold,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                vs::mag_sqr(grid.orientation[0]), 3.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
};

template <typename SrcTrait>
struct UpdateVelOp<::amr_wind_tests::FlatPlate, SrcTrait>
{
    void operator()(::amr_wind_tests::FlatPlate::DataType& data)
    {
        const auto& grid = data.grid();
        const auto& pos = grid.vel_pos;
        const auto& vel = grid.vel;
        amrex::Real rerr = 0.0_rt;
        for (int i = 0; i < grid.vel_pos.size(); ++i) {
            const amrex::Real val = pos[i].x() + pos[i].y() + pos[i].z();
            const vs::Vector vgold{val, val, val};
            rerr += vs::mag_sqr(vel[i] - vgold);
        }
        EXPECT_NEAR(
            rerr, 0.0_rt,
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }
};

template <typename SrcTrait>
struct ComputeForceOp<::amr_wind_tests::FlatPlate, SrcTrait>
{
    void operator()(::amr_wind_tests::FlatPlate::DataType& data)
    {
        auto& grid = data.grid();
        const int npts = data.meta().num_pts;
        const auto& dx = data.meta().dx;
        for (int ip = 0; ip < npts; ++ip) {
            const auto& tmat = grid.orientation[ip];
            // Effective velocity at the wing control point in local frame
            auto wvel = tmat & grid.vel[ip];
            // Set spanwise component to zero to get a pure 2D velocity
            wvel.y() = 0.0_rt;

            const auto vmag = vs::mag(wvel);
            const auto aoa = std::atan2(wvel.z(), wvel.x());

            // Make up some Cl, Cd values
            const auto cl = amr_wind::utils::two_pi() * aoa;
            const auto cd = cl * std::sin(aoa);

            // Assume unit chord
            const auto qval = 0.5_rt * vmag * vmag * dx[ip];
            const auto lift = qval * cl;
            const auto drag = qval * cd;
            // Determine unit vector parallel and perpendicular to velocity
            // vector
            const auto drag_dir = wvel.unit() & tmat;
            const auto lift_dir = drag_dir ^ tmat.y();

            // Compute force on fluid from this section of wing
            grid.force[ip] = -(lift_dir * lift + drag * drag_dir);
        }
    }
};

} // namespace ops

template class ::amr_wind::actuator::
    ActModel<::amr_wind_tests::FlatPlate, ::amr_wind::actuator::ActSrcLine>;
} // namespace amr_wind::actuator

namespace amr_wind_tests {

class ActPhysicsTest : public ::amr_wind::actuator::Actuator
{
public:
    explicit ActPhysicsTest(::amr_wind::CFDSim& sim)
        : ::amr_wind::actuator::Actuator(sim)
    {}

protected:
    void prepare_outputs() override {}
};

TEST_F(ActFlatPlateTest, act_model_init)
{
    initialize_mesh();
    sim().repo().declare_field("actuator_src_term", 3, 0);
    amr_wind::actuator::ActuatorContainer::ParticleType::NextID(1U);
    {
        amrex::ParmParse pp("Actuator.TestFlatPlateLine");
        pp.add("num_points", 11);
        pp.addarr(
            "start", amrex::Vector<amrex::Real>{16.0_rt, 12.0_rt, 16.0_rt});
        pp.addarr("end", amrex::Vector<amrex::Real>{16.0_rt, 20.0_rt, 16.0_rt});
        pp.addarr(
            "epsilon", amrex::Vector<amrex::Real>{3.0_rt, 3.0_rt, 3.0_rt});
        pp.add("pitch", 6.0_rt);
    }

    ::amr_wind::actuator::ActModel<FlatPlate> flat_plate(sim(), "F1", 0);
    {
        amr_wind::actuator::utils::ActParser pp(
            "Actuator.TestFlatPlateLine", "Actuator.F1");
        flat_plate.read_inputs(pp);
    }

    amrex::Vector<int> act_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    flat_plate.determine_root_proc(act_proc_count);
    flat_plate.init_actuator_source();

    const auto& info = flat_plate.info();
    EXPECT_EQ(info.root_proc, 0);
    EXPECT_EQ(info.procs.size(), amrex::ParallelDescriptor::NProcs());
}

TEST_F(ActFlatPlateTest, actuator_init)
{
    initialize_mesh();
    auto& vel = sim().repo().declare_field("velocity", 3, 3);
    auto& density = sim().repo().declare_field("density", 1, 3);
    density.setVal(1.0_rt);
    init_field(vel);
    amr_wind::actuator::ActuatorContainer::ParticleType::NextID(1U);

    amrex::Vector<std::string> actuators{"T1", "T2"};
    {
        amrex::ParmParse pp("Actuator");
        pp.addarr("labels", actuators);
        pp.add("type", std::string("TestFlatPlateLine"));
    }
    {
        amrex::ParmParse pp("Actuator.TestFlatPlateLine");
        pp.add("num_points", 11);
        pp.addarr(
            "epsilon", amrex::Vector<amrex::Real>{1.0_rt, 1.0_rt, 1.0_rt});
        pp.add("pitch", 6.0_rt);
    }
    {
        for (int i = 0; i < 2; ++i) {
            amrex::Real zloc = 8.0_rt + (16.0_rt * i);
            amrex::ParmParse pp("Actuator." + actuators[i]);
            pp.addarr(
                "start", amrex::Vector<amrex::Real>{16.0_rt, 12.0_rt, zloc});
            pp.addarr(
                "end", amrex::Vector<amrex::Real>{16.0_rt, 20.0_rt, zloc});
        }
    }

    ActPhysicsTest act(sim());
    act.pre_init_actions();
    act.post_init_actions();
}

TEST_F(ActFlatPlateTest, flat_plate_init)
{
    initialize_mesh();
    auto& vel = sim().repo().declare_field("velocity", 3, 3);
    auto& density = sim().repo().declare_field("density", 1, 3);
    density.setVal(1.0_rt);
    vel.setVal(10.0_rt, 0, 1, 3);
    amr_wind::actuator::ActuatorContainer::ParticleType::NextID(1U);

    amrex::Vector<std::string> actuators{"T1", "T2"};
    {
        amrex::ParmParse pp("Actuator");
        pp.addarr("labels", actuators);
        pp.add("type", std::string("FlatPlateLine"));
    }
    {
        amrex::ParmParse pp("Actuator.FlatPlateLine");
        pp.add("num_points", 11);
        pp.addarr(
            "epsilon", amrex::Vector<amrex::Real>{1.0_rt, 1.0_rt, 1.0_rt});
        pp.add("pitch", 6.0_rt);
    }
    {
        for (int i = 0; i < 2; ++i) {
            amrex::Real zloc = 8.0_rt + (16.0_rt * i);
            amrex::ParmParse pp("Actuator." + actuators[i]);
            pp.addarr(
                "start", amrex::Vector<amrex::Real>{16.0_rt, 12.0_rt, zloc});
            pp.addarr(
                "end", amrex::Vector<amrex::Real>{16.0_rt, 20.0_rt, zloc});
        }
    }

    ActPhysicsTest act(sim());
    act.pre_init_actions();
    act.post_init_actions();
}
} // namespace amr_wind_tests
