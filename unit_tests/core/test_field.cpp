#include "aw_test_utils/MeshTest.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class FieldRepoTest : public MeshTest
{
public:
    void declare_default_fields()
    {
        auto& frepo = mesh().field_repo();
        frepo.declare_field("vel", 3, 1, 2);
        frepo.declare_cc_field("density", 1);
        frepo.declare_nd_field("pressure", 1);
    }
};

TEST_F(FieldRepoTest, field_pre_declare)
{
    populate_parameters();
    create_mesh_instance();
    declare_default_fields();
    initialize_mesh();

    const auto& frepo = mesh().field_repo();
    EXPECT_EQ(frepo.num_fields(), 4);

    int num_cc = 0;
    int num_nd = 0;
    int num_invalid = 0;
    for (const auto& field : frepo.fields()) {
        switch (field->field_location()) {
        case amr_wind::FieldLoc::CELL:
            ++num_cc;
            break;

        case amr_wind::FieldLoc::NODE:
            ++num_nd;
            break;

        default:
            ++num_invalid;
            break;
        }
    }

    EXPECT_EQ(num_cc, 3);
    EXPECT_EQ(num_nd, 1);
    EXPECT_EQ(num_invalid, 0);
}

TEST_F(FieldRepoTest, field_declare_checks)
{
    initialize_mesh();

    auto& frepo = mesh().field_repo();
#if !(defined(AMREX_USE_MPI) && defined(__APPLE__))
    // Ensure that reserved field names aren't used
    EXPECT_THROW(
        frepo.declare_field("new_field__FS_Old", 1), amrex::RuntimeError);
    // Ensure check on allowable number of states
    EXPECT_THROW(
        frepo.declare_field("new_field", 1, 0, 5), amrex::RuntimeError);
#endif

    // Check that redeclaration returns the same field
    auto& velocity = frepo.declare_field("vel", 3);
    auto& vel1 = frepo.declare_field("vel", 3);
    EXPECT_EQ(&velocity, &vel1);

#if !(defined(AMREX_USE_MPI) && defined(__APPLE__))
    // Ensure that attempt to reregister field checks consistency
    EXPECT_THROW(frepo.declare_field("vel", 1, 1), amrex::RuntimeError);
#endif
}

TEST_F(FieldRepoTest, field_post_declare)
{
    initialize_mesh();
    declare_default_fields();

    const auto& frepo = mesh().field_repo();
    EXPECT_EQ(frepo.num_fields(), 4);

    // Check that the IDs are consistent
    {
        const auto& f1 = frepo.get_field("vel");
        const auto& f2 = frepo.get_field(0);
        EXPECT_EQ(f1.id(), f2.id());
    }
    // Velocity has two states...
    {
        const auto& f1 = frepo.get_field("density");
        const auto& f2 = frepo.get_field(2);
        EXPECT_EQ(f1.id(), f2.id());
    }
    {
        const auto& f1 = frepo.get_field("pressure");
        const auto& f2 = frepo.get_field(3);
        EXPECT_EQ(f1.id(), f2.id());
    }
}

TEST_F(FieldRepoTest, field_get)
{
    initialize_mesh();
    declare_default_fields();

    const auto& frepo = mesh().field_repo();
    auto& velf = frepo.get_field("vel");
    const auto& presf = frepo.get_field("pressure");
    const auto& vel_old = frepo.get_field("vel", amr_wind::FieldState::Old);

    EXPECT_EQ(velf.field_location(), amr_wind::FieldLoc::CELL);
    EXPECT_EQ(presf.field_location(), amr_wind::FieldLoc::NODE);
    EXPECT_EQ(vel_old.field_state(), amr_wind::FieldState::Old);

    const amrex::Real vx = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vy = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vz = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    velf.setVal(vx, 0);
    velf.setVal(vy, 1);
    velf.setVal(vz, 2);

    auto mf_vel = velf.vec_ptrs();
    amrex::Vector<amrex::Real> golds{{vx, vy, vz}};
    for (auto* it : mf_vel) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_vel = it->min(i);
            const auto max_vel = it->max(i);
            EXPECT_NEAR(
                min_vel, golds[i],
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                min_vel, max_vel,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
}

TEST_F(FieldRepoTest, field_multiple_states)
{
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& velocity = frepo.declare_field("vel", 3, 0, 2);
    auto& vel_old = velocity.state(amr_wind::FieldState::Old);
    auto& veldiff = frepo.declare_field("vel_diff", 3, 0);

    const amrex::Real vx = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vy = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vz = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));

    const amrex::Real vx_old = 20.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vy_old = 20.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vz_old = 20.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));

    velocity.setVal(amrex::Vector<amrex::Real>{vx, vy, vz});
    vel_old.setVal(amrex::Vector<amrex::Real>{vx_old, vy_old, vz_old});

    amr_wind::field_ops::lincomb(
        veldiff, 1.0_rt, velocity, 0, -1.0_rt, vel_old, 0, 0, 3, 0);

    const int nlevels = mesh().finestLevel() + 1;
    amrex::Vector<amrex::Real> golds{
        {(vx - vx_old), (vy - vy_old), (vz - vz_old)}};
    for (int lev = 0; lev < nlevels; ++lev) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = veldiff(lev).min(i);
            const auto max_val = veldiff(lev).max(i);
            EXPECT_NEAR(
                min_val, golds[i],
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                min_val, max_val,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
}

TEST_F(FieldRepoTest, field_location)
{
    initialize_mesh();

    auto& field_repo = mesh().field_repo();
    const auto& velocity = field_repo.declare_field("vel", 3, 0, 2);
    const auto& pressure =
        field_repo.declare_field("p", 1, 0, 1, amr_wind::FieldLoc::NODE);

    const auto& umac =
        field_repo.declare_field("umac", 1, 0, 1, amr_wind::FieldLoc::XFACE);
    const auto& vmac =
        field_repo.declare_field("vmac", 1, 0, 1, amr_wind::FieldLoc::YFACE);
    const auto& wmac =
        field_repo.declare_field("wmac", 1, 0, 1, amr_wind::FieldLoc::ZFACE);

    EXPECT_EQ(velocity.field_location(), amr_wind::FieldLoc::CELL);
    EXPECT_EQ(pressure.field_location(), amr_wind::FieldLoc::NODE);
    EXPECT_EQ(umac.field_location(), amr_wind::FieldLoc::XFACE);
    EXPECT_EQ(vmac.field_location(), amr_wind::FieldLoc::YFACE);
    EXPECT_EQ(wmac.field_location(), amr_wind::FieldLoc::ZFACE);

    const auto xf = amrex::IndexType(amrex::IntVect::TheDimensionVector(0));
    const auto yf = amrex::IndexType(amrex::IntVect::TheDimensionVector(1));
    const auto zf = amrex::IndexType(amrex::IntVect::TheDimensionVector(2));
    const int nlevels = field_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        EXPECT_EQ(velocity(lev).ixType(), amrex::IndexType::TheCellType());
        EXPECT_EQ(pressure(lev).ixType(), amrex::IndexType::TheNodeType());

        EXPECT_EQ(umac(lev).ixType(), xf);
        EXPECT_EQ(vmac(lev).ixType(), yf);
        EXPECT_EQ(wmac(lev).ixType(), zf);
    }
}

TEST_F(FieldRepoTest, field_advance_states)
{
    initialize_mesh();

    auto& field_repo = mesh().field_repo();
    auto& velocity = field_repo.declare_field("vel", 3, 0, 2);
    auto& vel_old = velocity.state(amr_wind::FieldState::Old);
    EXPECT_EQ(velocity.num_states(), 2);

    const amrex::Real vx = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vy = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));
    const amrex::Real vz = 30.0_rt + (5.0_rt * (amrex::Random() - 0.5_rt));

    velocity.setVal(amrex::Vector<amrex::Real>{vx, vy, vz});
    vel_old.setVal(std::numeric_limits<amrex::Real>::max());
    field_repo.advance_states();

    const int nlevels = field_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto old_min = vel_old(lev).min(i);
            const auto old_max = vel_old(lev).max(i);
            const auto new_min = velocity(lev).min(i);
            const auto new_max = velocity(lev).max(i);
            EXPECT_NEAR(
                old_min, new_min,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                old_max, new_max,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
}

TEST_F(FieldRepoTest, field_create_state)
{
    initialize_mesh();
    auto& field_repo = mesh().field_repo();
    auto& density = field_repo.declare_field("density");

    // Assert that there is no half state
    EXPECT_TRUE(density.num_states() == 1);
    EXPECT_FALSE(density.query_state(amr_wind::FieldState::NPH));

    // Create a new state and ensure that Field information is updated to
    // reflect this new state
    auto& rho_nph = density.create_state(amr_wind::FieldState::NPH);
    auto& rho_nph1 = density.state(amr_wind::FieldState::NPH);
    // Exact time states is still one
    EXPECT_TRUE(density.num_time_states() == 1);
    // But total states are two
    EXPECT_TRUE(density.num_states() == 2);
    EXPECT_TRUE(density.query_state(amr_wind::FieldState::NPH));
    EXPECT_EQ(&rho_nph, &rho_nph1);

    // Ensure that we can still "declare_field" and get the old field back as
    // long as the number of states is less than or equal to what has been
    // already created
    auto& rho_tmp = field_repo.declare_field("density");
    EXPECT_EQ(&rho_tmp, &density);

    const amrex::Real rho_val = 1.0_rt + (amrex::Random() - 0.5_rt);
    density.setVal(rho_val);
    amr_wind::field_ops::copy(
        rho_nph, density, 0, 0, density.num_comp(), density.num_grow());

    const int nlevels = field_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        EXPECT_NEAR(
            rho_nph(lev).min(0), rho_val,
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        EXPECT_NEAR(
            rho_nph(lev).max(0), rho_val,
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }

    // Ensure that advance states is a no op
    EXPECT_FALSE(density.query_state(amr_wind::FieldState::Old));
    density.advance_states();
}

TEST_F(FieldRepoTest, scratch_fields)
{
    populate_parameters();
    create_mesh_instance();

    auto& frepo = mesh().field_repo();
    auto& rho = frepo.declare_field("rho", 1, 0, 1);

    // Check that scratch field creation is disallowed before mesh is created
#if !(defined(AMREX_USE_MPI) && defined(__APPLE__))
    EXPECT_THROW(frepo.create_scratch_field(3, 0), amrex::RuntimeError);
#endif

    initialize_mesh();
    auto umac = frepo.create_scratch_field(3, 0, amr_wind::FieldLoc::XFACE);
    auto rho_nph = frepo.create_scratch_field(1, 0);

    const auto xf = amrex::IndexType(amrex::IntVect::TheDimensionVector(0));
    const int nlevels = frepo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        EXPECT_EQ((*umac)(lev).ixType(), xf);

        rho(lev).setVal(1.225_rt);
    }

    amr_wind::field_ops::copy(*rho_nph, rho, 0, 0, 1, 0);

    for (int lev = 0; lev < nlevels; ++lev) {
        EXPECT_NEAR(
            (*rho_nph)(lev).min(0), 1.225_rt,
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        EXPECT_NEAR(
            (*rho_nph)(lev).max(0), 1.225_rt,
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }
}

TEST_F(FieldRepoTest, int_scratch_fields)
{

    populate_parameters();
    create_mesh_instance();

    const auto& frepo = mesh().field_repo();
// Check that int scratch field creation is disallowed before mesh is created
#if !(defined(AMREX_USE_MPI) && defined(__APPLE__))
    EXPECT_THROW(
        frepo.create_int_scratch_field_on_host(1, 0), amrex::RuntimeError);
#endif
    initialize_mesh();

    auto ibcell_host = frepo.create_int_scratch_field_on_host(
        "iblank_cell_host", 1, 0, amr_wind::FieldLoc::CELL);
    auto ibnode_host = frepo.create_int_scratch_field_on_host(
        "iblank_node_host", 1, 0, amr_wind::FieldLoc::NODE);

    const int nlevels = frepo.num_active_levels();

    (*ibcell_host).setVal(1);
    (*ibnode_host).setVal(0);

    for (int lev = 0; lev < nlevels; ++lev) {
        EXPECT_EQ((*ibcell_host)(lev).max(0), 1);
        EXPECT_EQ((*ibcell_host)(lev).min(0), 1);
        EXPECT_EQ((*ibnode_host)(lev).max(0), 0);
        EXPECT_EQ((*ibnode_host)(lev).min(0), 0);
    }
}

TEST_F(FieldRepoTest, int_fields)
{
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& ibcell = frepo.declare_int_field("iblank_cell", 1, 0, 1);
    auto& ibnode = frepo.declare_int_field(
        "iblank_node", 1, 0, 1, amr_wind::FieldLoc::NODE);

    const int nlevels = frepo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        ibcell(lev).setVal(1);
        ibnode(lev).setVal(0);
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        EXPECT_EQ(ibcell(lev).max(0), 1);
        EXPECT_EQ(ibcell(lev).min(0), 1);
        EXPECT_EQ(ibnode(lev).max(0), 0);
        EXPECT_EQ(ibnode(lev).min(0), 0);
    }
}

TEST_F(FieldRepoTest, default_fillpatch_op)
{
    initialize_mesh();
    auto& frepo = mesh().field_repo();

    auto& velocity = frepo.declare_field("vel", 3, 1);
    EXPECT_FALSE(velocity.bc_initialized());
    EXPECT_FALSE(velocity.has_fillpatch_op());

    velocity.set_default_fillpatch_bc(sim().time());
    EXPECT_TRUE(velocity.bc_initialized());
    EXPECT_TRUE(velocity.has_fillpatch_op());

    velocity.setVal({10.0_rt, 20.0_rt, 30.0_rt}, 0);
    velocity.fillpatch(sim().time().current_time());
}

TEST_F(FieldRepoTest, field_subviews)
{
    initialize_mesh();
    auto& repo = mesh().field_repo();
    auto& velocity = repo.declare_field("vel", 3);

    velocity.setVal(10.0_rt, 0, 1);
    velocity.setVal(20.0_rt, 1, 1);
    velocity.setVal(30.0_rt, 2, 1);
    {
        auto wvel = velocity.subview(2);
        EXPECT_TRUE(wvel.num_comp() == 1);

        int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            EXPECT_NEAR(
                wvel(lev).min(0), 30.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                wvel(lev).max(0), 30.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }

    {
        auto vel2d = velocity.subview(0, 2);
        EXPECT_TRUE(vel2d.num_comp() == 2);
        int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            EXPECT_NEAR(
                vel2d(lev).min(0), 10.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                vel2d(lev).max(0), 10.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                vel2d(lev).min(1), 20.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
            EXPECT_NEAR(
                vel2d(lev).max(1), 20.0_rt,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
}

} // namespace amr_wind_tests
