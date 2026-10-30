#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/BCInterface.H"
#include "src/core/FieldFillPatchOps.H"
#include "AMReX_REAL.H"
#include "AMReX_Reduce.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

constexpr int nghost = 3;
constexpr amrex::Real temp_bc = 300.0_rt;

AMREX_GPU_HOST_DEVICE amrex::Real vel_interior(int i, int j, int k, int n)
{
    return 1.0_rt + (0.1_rt * i) + (0.01_rt * j) + (0.001_rt * k) +
           (0.25_rt * n);
}

AMREX_GPU_HOST_DEVICE amrex::Real temp_interior(int i, int j, int k)
{
    return 290.0_rt + (0.5_rt * i) + (0.05_rt * j) + (0.005_rt * k);
}

void init_interior(kynema_sgf::Field& velocity, kynema_sgf::Field& temperature)
{
    velocity.setVal(-99.0_rt, 0, velocity.num_comp(), nghost);
    temperature.setVal(-99.0_rt, 0, 1, nghost);
    for (int lev = 0; lev < velocity.repo().num_active_levels(); ++lev) {
        const auto& varrs = velocity(lev).arrays();
        const auto& tarrs = temperature(lev).arrays();
        amrex::ParallelFor(
            velocity(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    varrs[nbx](i, j, k, n) = vel_interior(i, j, k, n);
                }
                tarrs[nbx](i, j, k) = temp_interior(i, j, k);
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
    const bool check_velocity = true)
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
                                outflow ? vel_interior(iint, j, k, n)
                                        : ((n == 0) ? ux_bc : 0.0_rt);
                            err = amrex::max(
                                err, amrex::Math::abs(varr(i, j, k, n) - ref));
                        }
                    }
                    const amrex::Real tref =
                        outflow ? temp_interior(iint, j, k) : temp_bc;
                    return {amrex::max(
                        err, amrex::Math::abs(tarr(i, j, k) - tref))};
                });
        }
    }
    amrex::Real err = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
}

} // namespace

class MassInflowOutflowTest : public MeshTest
{
protected:
    void setup(
        const amrex::Real ux_bc,
        const bool high_velocity_pressure_outflow = false)
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{0, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        for (const std::string face : {"xlo", "xhi"}) {
            amrex::ParmParse pp(face);
            const bool high_face = face == "xhi";
            pp.add(
                "type", high_face && high_velocity_pressure_outflow
                            ? std::string("pressure_outflow")
                            : std::string("mass_inflow_outflow"));
            if (high_face && high_velocity_pressure_outflow) {
                pp.add("temperature_type", std::string("mass_inflow_outflow"));
            }
            pp.addarr("velocity", amrex::Vector<amrex::Real>{ux_bc, 0.0, 0.0});
            pp.addarr("temperature", amrex::Vector<amrex::Real>{temp_bc});
        }
        initialize_mesh();

        auto& repo = mesh().field_repo();
        m_vel = &repo.declare_field("velocity", AMREX_SPACEDIM, nghost, 2);
        m_temp = &repo.declare_field("temperature", 1, nghost, 1);

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
    EXPECT_EQ(ghost_error(*m_vel, *m_temp, ux_bc, false), 0.0_rt);
}

TEST_F(MassInflowOutflowTest, fillpatch_outflow_low)
{
    constexpr amrex::Real ux_bc = -1.0_rt;
    setup(ux_bc);
    m_vel->fillpatch(0.0_rt);
    m_temp->fillpatch(0.0_rt);
    EXPECT_EQ(ghost_error(*m_vel, *m_temp, ux_bc, true), 0.0_rt);
}

TEST_F(MassInflowOutflowTest, fillphysbc_outflow_low)
{
    constexpr amrex::Real ux_bc = -1.0_rt;
    setup(ux_bc);
    m_vel->fillphysbc(0.0_rt);
    m_temp->fillphysbc(0.0_rt);
    EXPECT_EQ(ghost_error(*m_vel, *m_temp, ux_bc, true), 0.0_rt);
}

TEST_F(MassInflowOutflowTest, scalar_outflow_uses_extrapolated_velocity)
{
    constexpr amrex::Real ux_bc = 1.0_rt;
    setup(ux_bc, true);
    m_vel->fillpatch(0.0_rt);
    m_temp->fillpatch(0.0_rt);

    EXPECT_EQ(ghost_error(*m_vel, *m_temp, ux_bc, false, false), 0.0_rt);
}

TEST_F(MassInflowOutflowTest, scalar_fills_preserve_velocity)
{
    constexpr amrex::Real ux_bc = 1.0_rt;
    setup(ux_bc);
    auto& velocity = (*m_vel)(0);
    amrex::MultiFab difference(
        velocity.boxArray(), velocity.DistributionMap(), velocity.nComp(),
        velocity.nGrowVect());

    for (const bool fillpatch : {true, false}) {
        SCOPED_TRACE(fillpatch ? "fillpatch" : "fillphysbc");
        amrex::MultiFab::Copy(
            difference, velocity, 0, 0, velocity.nComp(), velocity.nGrowVect());
        if (fillpatch) {
            m_temp->fillpatch(0.0_rt);
        } else {
            m_temp->fillphysbc(0.0_rt);
        }
        EXPECT_EQ(ghost_error(*m_vel, *m_temp, ux_bc, false, false), 0.0_rt);
        amrex::MultiFab::Subtract(
            difference, velocity, 0, 0, velocity.nComp(), velocity.nGrowVect());
        EXPECT_EQ(
            difference.norm0(0, velocity.nComp(), velocity.nGrowVect()),
            0.0_rt);
    }
}

} // namespace kynema_sgf_tests
