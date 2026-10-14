#include <numbers>

#include "ks_test_utils/MeshTest.H"

#include "src/wind_energy/actuator/FLLC.H"
#include "src/wind_energy/actuator/turbine/ActSrcDiskOp_Turbine.H"
#include "src/wind_energy/actuator/turbine/turbine_types.H"
#include "src/core/Slice.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

// Named (not anonymous) namespace: the device kernels of the source op are
// instantiated with these types, and nvcc requires external linkage there.
namespace turbine_disk {

namespace act = kynema_sgf::actuator;
namespace vs = kynema_sgf::vs;
namespace utils = kynema_sgf::utils;

//! Turbine actuator without an external solver, to drive the disk source op
struct DiskTurbine : public act::TurbineType
{
    using InfoType = act::TurbineInfo;
    using GridType = act::ActGrid;
    using MetaType = act::TurbineBaseData;
    using DataType = act::ActDataHolder<DiskTurbine>;

    static std::string identifier() { return "TestDiskTurbine"; }
};

using DiskSrcOp = act::ops::ActSrcOp<DiskTurbine, act::ActSrcDisk>;

constexpr int num_blades = 3;
constexpr int num_pts_blade = 8;
constexpr int num_pts_tower = 4;
// Radial spacing of the blade points (larger than 1 m so that dR and dR^2
// differ)
constexpr amrex::Real blade_dr = 3.0_rt;
const vs::Vector rotor_center{32.0_rt, 32.0_rt, 32.0_rt};

// Populate the actuator grid with the same layout as the external turbine
// models (hub, blades, tower) and set up the component views into it.
void init_disk_turbine(DiskTurbine::DataType& data)
{
    auto& info = data.info();
    auto& grid = data.grid();
    auto& meta = data.meta();

    info.bound_box =
        amrex::RealBox({0.0_rt, 0.0_rt, 0.0_rt}, {64.0_rt, 64.0_rt, 64.0_rt});

    meta.num_blades = num_blades;
    meta.num_pts_blade = num_pts_blade;
    meta.num_vel_pts_blade = num_pts_blade;
    meta.num_pts_tower = num_pts_tower;
    meta.rot_center = rotor_center;
    meta.rotor_frame = vs::Tensor::identity();

    const int npts = 1 + (num_blades * num_pts_blade) + num_pts_tower;
    grid.resize(npts);
    for (int ip = 0; ip < npts; ++ip) {
        grid.epsilon[ip] = vs::Vector(3.0_rt, 3.0_rt, 3.0_rt);
        grid.orientation[ip] = vs::Tensor::identity();
        grid.force[ip] = vs::Vector::zero();
    }

    grid.pos[0] = rotor_center;
    meta.hub.pos = utils::slice(grid.pos, 0, 1);
    meta.hub.force = utils::slice(grid.force, 0, 1);
    meta.hub.epsilon = utils::slice(grid.epsilon, 0, 1);
    meta.hub.orientation = utils::slice(grid.orientation, 0, 1);

    for (int ib = 0; ib < num_blades; ++ib) {
        const amrex::Real phi =
            2.0_rt * std::numbers::pi_v<amrex::Real> * ib / num_blades;
        const int start = 1 + (ib * num_pts_blade);
        for (int ip = 0; ip < num_pts_blade; ++ip) {
            const amrex::Real r = blade_dr * (ip + 1);
            grid.pos[start + ip] =
                rotor_center +
                vs::Vector(0.0_rt, r * std::cos(phi), r * std::sin(phi));
        }

        act::ComponentView cv;
        cv.pos = utils::slice(grid.pos, start, num_pts_blade);
        cv.force = utils::slice(grid.force, start, num_pts_blade);
        cv.epsilon = utils::slice(grid.epsilon, start, num_pts_blade);
        cv.orientation = utils::slice(grid.orientation, start, num_pts_blade);
        meta.blades.emplace_back(cv);
    }

    const int tstart = 1 + (num_blades * num_pts_blade);
    for (int ip = 0; ip < num_pts_tower; ++ip) {
        grid.pos[tstart + ip] =
            vs::Vector(36.0_rt, 32.0_rt, 16.0_rt + (3.0_rt * ip));
    }
    meta.tower.pos = utils::slice(grid.pos, tstart, num_pts_tower);
    meta.tower.force = utils::slice(grid.force, tstart, num_pts_tower);
    meta.tower.epsilon = utils::slice(grid.epsilon, tstart, num_pts_tower);
    meta.tower.orientation =
        utils::slice(grid.orientation, tstart, num_pts_tower);
}

// Spread the actuator forces with an already set up source op
void spread_source(
    DiskSrcOp& op, const amrex::Geometry& geom, kynema_sgf::Field& src)
{
    src.setVal(0.0_rt);
    for (amrex::MFIter mfi(src(0)); mfi.isValid(); ++mfi) {
        op(0, mfi, geom);
    }
    amrex::Gpu::streamSynchronize();
}

// Largest difference between two source fields, per component
amrex::Real
max_diff(const amrex::MultiFab& lhs, const amrex::MultiFab& rhs, const int comp)
{
    amrex::MultiFab diff(lhs.boxArray(), lhs.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(diff, lhs, comp, 0, 1, 0);
    amrex::MultiFab::Subtract(diff, rhs, comp, 0, 1, 0);
    return diff.norm0(0);
}

} // namespace turbine_disk

class TurbineDiskSrcTest : public MeshTest
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
            amrex::Vector<amrex::Real> probhi{{64.0_rt, 64.0_rt, 64.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

// The spreading kernel must read only the device copy of the actuator data
// made in setup_op. The host arrays (and the component views into them) are
// not accessible from a GPU kernel (kynema/kynema-sgf#1977). Overwriting the
// host data after setup_op must therefore leave the source unchanged; a
// kernel that dereferences host pointers picks up the new values on a CPU
// build and fails here, instead of passing silently as it does on CPU.
TEST_F(TurbineDiskSrcTest, spreads_device_copy_of_actuator_data)
{
    namespace td = turbine_disk;
    initialize_mesh();
    auto& src = sim().repo().declare_field("actuator_src_term", 3, 0);
    const auto& geom = sim().mesh().Geom(0);

    td::DiskTurbine::DataType data(sim(), "disk", 0);
    td::init_disk_turbine(data);

    auto& grid = data.grid();
    for (int ip = 0; ip < grid.force.size(); ++ip) {
        grid.force[ip] =
            td::vs::Vector(1.0_rt + (0.1_rt * ip), 0.2_rt, -0.1_rt * (ip % 3));
    }

    td::DiskSrcOp op(data);
    op.initialize();
    op.setup_op();
    td::spread_source(op, geom, src);

    amrex::MultiFab ref(src(0).boxArray(), src(0).DistributionMap(), 3, 0);
    amrex::MultiFab::Copy(ref, src(0), 0, 0, 3, 0);
    ASSERT_FALSE(ref.contains_nan());
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        ASSERT_GT(ref.norm0(n), 0.0_rt) << "component " << n;
    }

    // Change the host data without copying it to the device again
    for (int ip = 0; ip < grid.force.size(); ++ip) {
        grid.force[ip] = td::vs::Vector(-1.0e3_rt, 1.0e3_rt, 1.0e3_rt);
        grid.epsilon[ip] = td::vs::Vector(0.5_rt, 0.5_rt, 0.5_rt);
        grid.orientation[ip] = td::vs::Tensor::zero();
    }
    td::spread_source(op, geom, src);

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    EXPECT_FALSE(src(0).contains_nan());
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        EXPECT_LE(td::max_diff(src(0), ref, n), tol * ref.norm0(n))
            << "component " << n;
    }

    // Once copied to the device the new host data does change the source
    op.setup_op();
    td::spread_source(op, geom, src);
    EXPECT_FALSE(src(0).contains_nan());
    EXPECT_GT(td::max_diff(src(0), ref, 0), 0.1_rt * ref.norm0(0));
}

} // namespace kynema_sgf_tests
