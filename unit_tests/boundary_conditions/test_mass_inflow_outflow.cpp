#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/BCInterface.H"
#include "src/core/FieldFillPatchOps.H"
#include "src/incflo.H"
#include "AMReX_REAL.H"
#include "AMReX_Reduce.H"
#include <limits>
#include <sstream>

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

constexpr int nghost = 3;
constexpr amrex::Real temp_bc = 300.0_rt;
constexpr amrex::Real tolerance =
    std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

struct BoundaryFillStats
{
    int calls{0};
    amrex::Long cells{0};
    kynema_sgf::FieldState state{kynema_sgf::FieldState::New};
};

class StateBoundaryFill
    : public kynema_sgf::FieldFillPatchOps<kynema_sgf::FieldBCDirichlet>
{
public:
    StateBoundaryFill(
        kynema_sgf::Field& field,
        const amrex::AmrCore& mesh,
        const kynema_sgf::SimTime& time,
        BoundaryFillStats& stats)
        : FieldFillPatchOps(field, mesh, time), m_stats(stats)
    {}

    void fillphysbc(
        int lev,
        amrex::Real /*time*/,
        amrex::MultiFab& mfab,
        const amrex::IntVect& /*nghost*/,
        kynema_sgf::FieldState state) override
    {
        ++m_stats.calls;
        m_stats.state = state;
        for (int idx = 0; idx < mfab.boxArray().size(); ++idx) {
            EXPECT_EQ(mfab.boxArray()[idx].length(0), 1);
        }
        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            const int source_idx =
                kynema_sgf::field_impl::matching_grid_index(mfi, m_field(lev));
            EXPECT_TRUE(
                m_field(lev).boxArray()[source_idx].contains(mfi.validbox()));
            m_stats.cells += mfab[mfi].box().numPts();
        }
        mfab.setVal(state == kynema_sgf::FieldState::New ? 1.0_rt : -1.0_rt);
    }

private:
    BoundaryFillStats& m_stats;
};

AMREX_GPU_HOST_DEVICE amrex::Real vel_interior(int i, int j, int k, int n)
{
    return 1.0_rt + (0.1_rt * i) + (0.01_rt * j) + (0.001_rt * k) +
           (0.25_rt * n);
}

AMREX_GPU_HOST_DEVICE amrex::Real temp_interior(int i, int j, int k)
{
    return 290.0_rt + (0.5_rt * i) + (0.05_rt * j) + (0.005_rt * k);
}

void init_interior(
    kynema_sgf::Field& velocity,
    kynema_sgf::Field& temperature,
    const amrex::Real offset = 0.0_rt)
{
    velocity.setVal(-99.0_rt, 0, velocity.num_comp(), nghost);
    temperature.setVal(-99.0_rt, 0, 1, nghost);
    for (int lev = 0; lev < velocity.repo().num_active_levels(); ++lev) {
        const auto& varrs = velocity(lev).arrays();
        const auto& tarrs = temperature(lev).arrays();
        amrex::ParallelFor(
            velocity(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    varrs[nbx](i, j, k, n) = vel_interior(i, j, k, n) + offset;
                }
                tarrs[nbx](i, j, k) = temp_interior(i, j, k) + offset;
            });
    }
    amrex::Gpu::streamSynchronize();
}

/** Max error in the x-boundary ghost layers. On the outflow side the ghost
 *  values must match the adjacent interior cell, on the inflow side they must
 *  match the prescribed Dirichlet values.
 */
amrex::Real ghost_error(
    const kynema_sgf::Field& velocity,
    const kynema_sgf::Field& temperature,
    const amrex::Real ux_bc,
    const bool lo_outflow,
    const bool check_velocity = true,
    const amrex::Real offset = 0.0_rt,
    const amrex::Real temperature_bc = temp_bc)
{
    const auto& domain = velocity.repo().mesh().Geom(0).Domain();
    const int ilo = domain.smallEnd(0);
    const int ihi = domain.bigEnd(0);
    amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (amrex::MFIter mfi(velocity(0)); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto varr = velocity(0).const_array(mfi);
        const auto tarr = temperature(0).const_array(mfi);
        for (const bool is_low : {true, false}) {
            const bool outflow = (is_low == lo_outflow);
            const int iint = is_low ? ilo : ihi;
            if ((is_low && vbx.smallEnd(0) != ilo) ||
                (!is_low && vbx.bigEnd(0) != ihi)) {
                continue;
            }
            const auto bx = is_low ? amrex::adjCellLo(vbx, 0, nghost)
                                   : amrex::adjCellHi(vbx, 0, nghost);
            reduce_op.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(
                    int i, int j, int k) noexcept -> ReduceTuple {
                    amrex::Real err = 0.0_rt;
                    if (check_velocity) {
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            const amrex::Real ref =
                                outflow ? vel_interior(iint, j, k, n) + offset
                                        : ((n == 0) ? ux_bc : 0.0_rt);
                            const auto delta =
                                amrex::Math::abs(varr(i, j, k, n) - ref);
                            err = amrex::max(
                                err,
                                std::isfinite(delta)
                                    ? delta
                                    : std::numeric_limits<amrex::Real>::max());
                        }
                    }
                    const amrex::Real tref =
                        outflow ? temp_interior(iint, j, k) + offset
                                : temperature_bc;
                    const auto delta = amrex::Math::abs(tarr(i, j, k) - tref);
                    return {amrex::max(
                        err, std::isfinite(delta)
                                 ? delta
                                 : std::numeric_limits<amrex::Real>::max())};
                });
        }
    }
    amrex::Real err = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

amrex::Long boundary_ghost_errors(
    const kynema_sgf::Field& field,
    const int lev,
    const amrex::Real inflow_value,
    const bool low_outflow = false)
{
    const auto domain = field.repo().mesh().Geom(lev).Domain();
    const int ncomp = field.num_comp();
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<amrex::Long> data(op);
    using Tuple = typename decltype(data)::Type;
    for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
        const auto arr = field(lev).const_array(mfi);
        const auto valid = mfi.validbox();
        for (int side = 0; side < 2; ++side) {
            const int interior =
                side == 0 ? domain.smallEnd(0) : domain.bigEnd(0);
            if ((side == 0 ? valid.smallEnd(0) : valid.bigEnd(0)) != interior) {
                continue;
            }
            const auto box =
                side == 0 ? amrex::adjCellLo(valid, 0, field.num_grow()[0])
                          : amrex::adjCellHi(valid, 0, field.num_grow()[0]);
            op.eval(
                box, data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> Tuple {
                    amrex::Long errors = 0;
                    for (int n = 0; n < ncomp; ++n) {
                        const amrex::Real expected =
                            (side == 0) == low_outflow
                                ? arr(interior, j, k, n)
                                : (n == 0 ? inflow_value : 0.0_rt);
                        errors += !(
                            std::abs(arr(i, j, k, n) - expected) <= tolerance);
                    }
                    return {errors};
                });
        }
    }
    auto errors = amrex::get<0>(data.value());
    amrex::ParallelDescriptor::ReduceLongSum(errors);
    return errors;
}

} // namespace

class MassInflowOutflowTest : public MeshTest
{
protected:
    void setup(const amrex::Real ux_bc, const bool refined = false)
    {
        MeshTest::populate_parameters();
        amrex::ParmParse("amr").addarr(
            "n_cell", amrex::Vector<int>{16, 16, 16});
        amrex::ParmParse("amr").add("max_grid_size", 8);
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{0, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        for (const std::string face : {"xlo", "xhi"}) {
            amrex::ParmParse pp(face);
            pp.add("type", std::string("mass_inflow_outflow"));
            pp.addarr("velocity", amrex::Vector<amrex::Real>{ux_bc, 0.0, 0.0});
            pp.addarr("temperature", amrex::Vector<amrex::Real>{temp_bc});
        }
        if (refined) {
            amrex::ParmParse("amr").add("max_level", 1);
            amrex::ParmParse("amr").add("blocking_factor", 2);
            create_mesh_instance<RefineMesh>();
            auto refinement =
                std::make_unique<kynema_sgf::CartBoxRefinement>(sim());
            std::stringstream boxes("1\n2\n0 0 0 2 8 8\n6 0 0 8 8 8\n");
            refinement->read_inputs(mesh(), boxes);
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(refinement));
        }
        initialize_mesh();

        auto& repo = mesh().field_repo();
        m_vel = &repo.declare_field("velocity", AMREX_SPACEDIM, nghost, 2);
        m_temp = &repo.declare_field("temperature", 1, nghost, 2);

        m_vel->register_fill_patch_op<
            kynema_sgf::FieldFillPatchOps<kynema_sgf::FieldBCDirichlet>>(
            mesh(), time(), kynema_sgf::FieldInterpolator::CellConsLinear);
        m_temp->register_fill_patch_op<
            kynema_sgf::FieldFillPatchOps<kynema_sgf::FieldBCDirichlet>>(
            mesh(), time(), kynema_sgf::FieldInterpolator::CellConsLinear);
        kynema_sgf::BCVelocity bc_vel(*m_vel);
        bc_vel();
        kynema_sgf::BCScalar bc_temp(*m_temp);
        bc_temp();

        init_interior(*m_vel, *m_temp);
    }

    kynema_sgf::Field* m_vel{nullptr};
    kynema_sgf::Field* m_temp{nullptr};
};

TEST_F(MassInflowOutflowTest, bcrec_is_direction_dependent)
{
    setup(1.0_rt);
    EXPECT_TRUE(m_vel->has_inout_bndry());
    EXPECT_TRUE(m_vel->state(kynema_sgf::FieldState::Old).has_inout_bndry());
    for (const auto* fld : {m_vel, m_temp}) {
        for (const auto& bcr : fld->bcrec()) {
            EXPECT_EQ(bcr.lo(0), amrex::BCType::direction_dependent);
            EXPECT_EQ(bcr.hi(0), amrex::BCType::direction_dependent);
            EXPECT_EQ(bcr.lo(1), amrex::BCType::int_dir);
        }
    }
}

TEST_F(MassInflowOutflowTest, fillpatch_outflow_high)
{
    constexpr amrex::Real ux_bc = 1.0_rt;
    setup(ux_bc);
    m_vel->fillpatch(0.0_rt);
    m_temp->fillpatch(0.0_rt);
    EXPECT_NEAR(ghost_error(*m_vel, *m_temp, ux_bc, false), 0.0_rt, tolerance);
}

TEST_F(MassInflowOutflowTest, fillpatch_outflow_low)
{
    constexpr amrex::Real ux_bc = -1.0_rt;
    setup(ux_bc);
    m_vel->fillpatch(0.0_rt);
    m_temp->fillpatch(0.0_rt);
    EXPECT_NEAR(ghost_error(*m_vel, *m_temp, ux_bc, true), 0.0_rt, tolerance);
}

TEST_F(MassInflowOutflowTest, fillphysbc_outflow_low)
{
    constexpr amrex::Real ux_bc = -1.0_rt;
    setup(ux_bc);
    m_vel->fillphysbc(0.0_rt);
    m_temp->fillphysbc(0.0_rt);
    EXPECT_NEAR(ghost_error(*m_vel, *m_temp, ux_bc, true), 0.0_rt, tolerance);
}

TEST_F(MassInflowOutflowTest, scalar_uses_matching_state_and_boundary_strips)
{
    setup(1.0_rt);
    BoundaryFillStats stats;
    m_vel->register_fill_patch_op<StateBoundaryFill>(mesh(), time(), stats);
    auto& old_velocity = m_vel->state(kynema_sgf::FieldState::Old);
    auto& old_temperature = m_temp->state(kynema_sgf::FieldState::Old);
    init_interior(old_velocity, old_temperature);
    amrex::MultiFab difference(
        (*m_vel)(0).boxArray(), (*m_vel)(0).DistributionMap(), 3, nghost);
    amrex::MultiFab::Copy(difference, (*m_vel)(0), 0, 0, 3, nghost);
    old_temperature.fillphysbc(0.0_rt);
    EXPECT_EQ(stats.state, kynema_sgf::FieldState::Old);
    EXPECT_EQ(stats.calls, 2);
    EXPECT_NEAR(
        ghost_error(old_velocity, old_temperature, -1.0_rt, true, false),
        0.0_rt, tolerance);
    EXPECT_LT(
        stats.cells, amrex::grow(mesh().Geom(0).Domain(), nghost).numPts());
    amrex::MultiFab::Subtract(difference, (*m_vel)(0), 0, 0, 3, nghost);
    for (int n = 0; n < 3; ++n) {
        EXPECT_NEAR(difference.norm0(n, nghost), 0.0_rt, tolerance);
    }

    m_temp->fillpatch(0.0_rt);
    EXPECT_EQ(stats.state, kynema_sgf::FieldState::New);
    EXPECT_NEAR(
        ghost_error(*m_vel, *m_temp, 1.0_rt, false, false), 0.0_rt, tolerance);
}

TEST_F(MassInflowOutflowTest, interior_grids_do_not_evaluate_velocity)
{
    setup(1.0_rt);
    BoundaryFillStats stats;
    m_vel->register_fill_patch_op<StateBoundaryFill>(mesh(), time(), stats);
    auto interior = mesh().Geom(0).Domain();
    interior.grow(0, -4);
    ASSERT_TRUE(interior.ok());
    const amrex::BoxArray ba(interior);
    amrex::MultiFab scratch(ba, amrex::DistributionMapping(ba), 1, nghost);
    scratch.setVal(123.0_rt);
    m_temp->fill_inout_outflow(0, 0.0_rt, scratch, amrex::IntVect(nghost));
    EXPECT_EQ(stats.calls, 0);
    EXPECT_NEAR(scratch.norm0(0, nghost), 123.0_rt, tolerance);
    m_temp->fill_inout_outflow(0, 0.0_rt, (*m_temp)(0), amrex::IntVect(0));
    EXPECT_EQ(stats.calls, 0);
}

TEST_F(MassInflowOutflowTest, missing_velocity_state_uses_new)
{
    setup(1.0_rt);
    BoundaryFillStats stats;
    m_vel->register_fill_patch_op<StateBoundaryFill>(mesh(), time(), stats);
    auto& half_temperature = m_temp->create_state(kynema_sgf::FieldState::NPH);
    half_temperature.setVal(temp_bc);
    half_temperature.fillphysbc(0.0_rt);
    EXPECT_EQ(stats.state, kynema_sgf::FieldState::New);
}

TEST_F(MassInflowOutflowTest, scalar_uses_existing_half_time_velocity)
{
    setup(1.0_rt);
    BoundaryFillStats stats;
    m_vel->register_fill_patch_op<StateBoundaryFill>(mesh(), time(), stats);
    auto& half_velocity = m_vel->create_state(kynema_sgf::FieldState::NPH);
    auto& half_temperature = m_temp->create_state(kynema_sgf::FieldState::NPH);
    init_interior(half_velocity, half_temperature, 10.0_rt);
    amrex::MultiFab original(
        (*m_temp)(0).boxArray(), (*m_temp)(0).DistributionMap(), 1, nghost);
    amrex::MultiFab::Copy(original, (*m_temp)(0), 0, 0, 1, nghost);
    half_temperature.fillphysbc(0.05_rt);
    EXPECT_EQ(stats.state, kynema_sgf::FieldState::NPH);
    EXPECT_NEAR(
        ghost_error(
            half_velocity, half_temperature, -1.0_rt, true, false, 10.0_rt),
        0.0_rt, tolerance);
    amrex::MultiFab::Subtract(original, (*m_temp)(0), 0, 0, 1, nghost);
    EXPECT_NEAR(original.norm0(0, nghost), 0.0_rt, tolerance);
}

TEST_F(MassInflowOutflowTest, successive_fills_refresh_boundary_values)
{
    setup(1.0_rt);
    for (const bool physical_only : {false, true}) {
        SCOPED_TRACE(physical_only ? "fillphysbc" : "fillpatch");
        for (int step = 0; step < 2; ++step) {
            const amrex::Real ux = step == 0 ? 1.0_rt : -1.0_rt;
            const amrex::Real offset = 20.0_rt * step;
            const amrex::Real theta = temp_bc + offset;
            for (const int side : {0, 1}) {
                const amrex::Orientation ori(
                    0, side == 0 ? amrex::Orientation::low
                                 : amrex::Orientation::high);
                m_vel->bc_values()[ori][0] = ux;
                m_temp->bc_values()[ori][0] = theta;
            }
            m_vel->copy_bc_to_device();
            m_temp->copy_bc_to_device();
            init_interior(*m_vel, *m_temp, offset);
            if (physical_only) {
                m_vel->fillphysbc(0.1_rt * step);
                m_temp->fillphysbc(0.1_rt * step);
            } else {
                m_vel->fillpatch(0.1_rt * step);
                m_temp->fillpatch(0.1_rt * step);
            }
            EXPECT_NEAR(
                ghost_error(
                    *m_vel, *m_temp, ux, step == 1, true, offset, theta),
                0.0_rt, tolerance);
        }
    }
}

TEST_F(MassInflowOutflowTest, coarse_to_fine_fills_complete_outflow)
{
    setup(1.0_rt, true);
    ASSERT_EQ(mesh().num_levels(), 2);
    for (const amrex::Real ux : {1.0_rt, -1.0_rt}) {
        for (const int side : {0, 1}) {
            m_vel->bc_values()[amrex::Orientation(
                0, side == 0 ? amrex::Orientation::low
                             : amrex::Orientation::high)][0] = ux;
        }
        m_vel->copy_bc_to_device();
        m_vel->setVal(2.0_rt);
        m_temp->setVal(290.0_rt);
        (*m_vel)(1).setVal(-99.0_rt);
        (*m_temp)(1).setVal(-99.0_rt);
        m_vel->fillpatch_from_coarse(1, 0.0_rt, (*m_vel)(1), nghost);
        m_temp->fillpatch_from_coarse(1, 0.0_rt, (*m_temp)(1), nghost);
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_NEAR((*m_vel)(1).min(n), 2.0_rt, tolerance);
            EXPECT_NEAR((*m_vel)(1).max(n), 2.0_rt, tolerance);
        }
        EXPECT_NEAR((*m_temp)(1).min(0), 290.0_rt, tolerance);
        EXPECT_NEAR((*m_temp)(1).max(0), 290.0_rt, tolerance);
        EXPECT_EQ(boundary_ghost_errors(*m_vel, 1, ux, ux < 0.0_rt), 0);
        EXPECT_EQ(boundary_ghost_errors(*m_temp, 1, temp_bc, ux < 0.0_rt), 0);
    }
}

TEST_F(MassInflowOutflowTest, amr_timestep_refreshes_field_states)
{
    MeshTest::populate_parameters();
    amrex::ParmParse("amr").addarr("n_cell", amrex::Vector<int>{16, 16, 8});
    amrex::ParmParse("amr").add("n_error_buf", 0);
    amrex::ParmParse("amr").add("max_level", 1);
    amrex::ParmParse("amr").add("blocking_factor", 2);
    amrex::ParmParse("amr").add("max_grid_size", 8);
    amrex::ParmParse("geometry")
        .addarr("is_periodic", amrex::Vector<int>{0, 1, 1});
    for (const auto* face : {"xlo", "xhi"}) {
        amrex::ParmParse pp(face);
        pp.add("type", std::string("mass_inflow_outflow"));
        pp.addarr(
            "velocity", amrex::Vector<amrex::Real>{1.0_rt, 0.0_rt, 0.0_rt});
        pp.addarr("density", amrex::Vector<amrex::Real>{1.0_rt});
        pp.addarr("temperature", amrex::Vector<amrex::Real>{temp_bc});
    }
    amrex::ParmParse("incflo").addarr(
        "physics", amrex::Vector<std::string>{"FreeStream"});
    amrex::ParmParse("incflo").add("use_godunov", 1);
    amrex::ParmParse("incflo").add("do_initial_proj", false);
    amrex::ParmParse("incflo").add("initial_iterations", 0);
    amrex::ParmParse("transport").add("viscosity", 0.0_rt);
    amrex::ParmParse("turbulence").add("model", std::string("Laminar"));
    amrex::ParmParse("incflo").addarr(
        "velocity", amrex::Vector<amrex::Real>{1.0_rt, 0.0_rt, 0.0_rt});
    amrex::ParmParse("incflo").add("density", 1.0_rt);
    amrex::ParmParse("incflo").add("temperature", 290.0_rt);
    amrex::ParmParse("tagging").addarr(
        "labels", amrex::Vector<std::string>{"refined"});
    amrex::ParmParse pp("tagging.refined");
    pp.add("type", std::string("GeometryRefinement"));
    pp.addarr("shapes", amrex::Vector<std::string>{"box"});
    pp.add("level", 0);
    amrex::ParmParse box("tagging.refined.box");
    box.add("type", std::string("box"));
    box.addarr("origin", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
    box.addarr("xaxis", amrex::Vector<amrex::Real>{8.0_rt, 0.0_rt, 0.0_rt});
    box.addarr("yaxis", amrex::Vector<amrex::Real>{0.0_rt, 4.0_rt, 0.0_rt});
    box.addarr("zaxis", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 8.0_rt});
    amrex::ParmParse("time").add("plot_interval", -1);
    amrex::ParmParse("time").add("checkpoint_interval", -1);
    amrex::ParmParse("time").add("regrid_interval", -1);
    reset_prob_domain();
    incflo driver;
    driver.sim().pde_manager().register_transport_pde("Temperature");
    auto& velocity = driver.velocity();
    auto& temperature = driver.temperature();
    velocity.create_state(kynema_sgf::FieldState::NPH);
    temperature.create_state(kynema_sgf::FieldState::NPH);
    driver.InitData();
    ASSERT_EQ(driver.repo().num_active_levels(), 2);
    ASSERT_LT(driver.boxArray(1).numPts(), driver.Geom(1).Domain().numPts());
    for (int step = 0; step < 2; ++step) {
        SCOPED_TRACE(step);
        velocity.state(kynema_sgf::FieldState::NPH).setVal(-99.0_rt);
        temperature.state(kynema_sgf::FieldState::NPH).setVal(-99.0_rt);
        driver.compute_dt();
        driver.pre_advance_stage1();
        driver.pre_advance_stage2();
        driver.prepare_time_step();
        for (int lev = 0; lev < 2; ++lev) {
            for (auto* field : {&velocity, &temperature}) {
                amrex::MultiFab difference(
                    (*field)(lev).boxArray(), (*field)(lev).DistributionMap(),
                    field->num_comp(), 0);
                amrex::MultiFab::Copy(
                    difference, field->state(kynema_sgf::FieldState::NPH)(lev),
                    0, 0, field->num_comp(), 0);
                amrex::MultiFab::Subtract(
                    difference, (*field)(lev), 0, 0, field->num_comp(), 0);
                for (int n = 0; n < field->num_comp(); ++n) {
                    EXPECT_NEAR(difference.norm0(n), 0.0_rt, tolerance);
                }
            }
            EXPECT_EQ(
                boundary_ghost_errors(
                    velocity.state(kynema_sgf::FieldState::NPH), lev, 1.0_rt),
                0);
            EXPECT_EQ(
                boundary_ghost_errors(
                    temperature.state(kynema_sgf::FieldState::NPH), lev,
                    temp_bc),
                0);
        }
        driver.ApplyPredictor();
        driver.post_advance_work();
        for (int lev = 0; lev < 2; ++lev) {
            for (const auto state :
                 {kynema_sgf::FieldState::Old, kynema_sgf::FieldState::New}) {
                EXPECT_EQ(
                    boundary_ghost_errors(velocity.state(state), lev, 1.0_rt),
                    0);
                EXPECT_EQ(
                    boundary_ghost_errors(
                        temperature.state(state), lev, temp_bc),
                    0);
            }
        }
        EXPECT_GT(temperature(0).max(0), 290.0_rt + tolerance);
    }
}

TEST_F(MassInflowOutflowTest, ghosts_reaching_boundary_without_valid_face)
{
    setup(-1.0_rt);
    BoundaryFillStats stats;
    m_vel->register_fill_patch_op<StateBoundaryFill>(mesh(), time(), stats);
    auto box = mesh().boxArray(0)[0];
    box.setRange(0, 2, 2);
    const amrex::BoxArray ba(box);
    amrex::MultiFab scratch(ba, amrex::DistributionMapping(ba), 1, nghost);
    scratch.setVal(123.0_rt);
    scratch.setVal(
        -99.0_rt, amrex::adjCellLo(mesh().Geom(0).Domain(), 0, 1), 0, 1,
        nghost);
    m_temp->state(kynema_sgf::FieldState::Old)
        .fill_inout_outflow(0, 0.0_rt, scratch, amrex::IntVect(nghost));
    EXPECT_EQ(stats.calls, 1);
    const auto ghosts = amrex::grow(box, nghost) &
                        amrex::adjCellLo(mesh().Geom(0).Domain(), 0, 1);
    EXPECT_NEAR(scratch.min(ghosts, 0, nghost), 123.0_rt, tolerance);
}

} // namespace kynema_sgf_tests
