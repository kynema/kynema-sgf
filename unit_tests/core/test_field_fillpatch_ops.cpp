#include "aw_test_utils/MeshTest.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/projection/nodal_projection_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

struct TestProfile
{
    struct DeviceOp
    {
        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::Array4<amrex::Real> const& field,
            amrex::GeometryData const& /*unused*/,
            const amrex::Real /*unused*/,
            amrex::Orientation /*unused*/,
            const int comp,
            const int dcomp,
            const int orig_comp) const
        {
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> vel{
                1.0_rt, 2.0_rt, 3.0_rt};
            field(iv[0], iv[1], iv[2], dcomp + comp) = vel[orig_comp + comp];
        }
    };

    using DeviceType = DeviceOp;

    static std::string identifier() { return "TestProfile"; }

    explicit TestProfile(const amr_wind::Field& /*unused*/) {}

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

amrex::Real get_field_err(
    amr_wind::Field& field, const bool check_all_ghost, const int comp = 0)
{
    const int lev = 0;
    amrex::Real error_total = 0;
    const int ncomp = field.num_comp();

    error_total += amrex::ReduceSum(
        field(lev), field(lev).nGrow(),
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) {
                // Do indices manually to check against box operations in code
                if (ncomp == 1) {
                    // Checking a MAC velocity
                    if (comp == 0) {
                        // Checking u
                        if (check_all_ghost) {
                            // Ghost cells in boundaries
                            if (i == -1 || j == -1 || k == -1 || i == 9 ||
                                j == 8 || k == 8) {
                                error += std::abs(f_arr(i, j, k) - 1.0_rt);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0_rt);
                            }
                        } else {
                            // Valid face-normal cells
                            if ((i == 0 || i == 8) && (j >= 0 && j <= 7) &&
                                (k >= 0 && k <= 7)) {
                                error += std::abs(f_arr(i, j, k) - 1.0_rt);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0_rt);
                            }
                        }
                    } else if (comp == 1) {
                        // Checking v
                        if (check_all_ghost) {
                            // Ghost cells in boundaries
                            if (i == -1 || j == -1 || k == -1 || i == 8 ||
                                j == 9 || k == 8) {
                                error += std::abs(f_arr(i, j, k) - 2.0_rt);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0_rt);
                            }
                        } else {
                            // Valid face-normal cells
                            if ((j == 0 || j == 8) && (i >= 0 && i <= 7) &&
                                (k >= 0 && k <= 7)) {
                                error += std::abs(f_arr(i, j, k) - 2.0_rt);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0_rt);
                            }
                        }
                    } else {
                        // Checking w
                        if (check_all_ghost) {
                            // Ghost cells in boundaries
                            if (i == -1 || j == -1 || k == -1 || i == 8 ||
                                j == 8 || k == 9) {
                                error += std::abs(f_arr(i, j, k) - 3.0_rt);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0_rt);
                            }
                        } else {
                            // Valid face-normal cells
                            if ((k == 0 || k == 8) && (i >= 0 && i <= 7) &&
                                (j >= 0 && j <= 7)) {
                                error += std::abs(f_arr(i, j, k) - 3.0_rt);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0_rt);
                            }
                        }
                    }
                } else {
                    // Checking cell-centered velocity
                    if (check_all_ghost) {
                        // Ghost cells in boundaries
                        if (i <= -1 || j <= -1 || k <= -1 || i >= 8 || j >= 8 ||
                            k >= 8) {
                            error += std::abs(f_arr(i, j, k, 0) - 1.0_rt);
                            error += std::abs(f_arr(i, j, k, 1) - 2.0_rt);
                            error += std::abs(f_arr(i, j, k, 2) - 3.0_rt);
                        } else {
                            error += std::abs(f_arr(i, j, k, 0) - 0.0_rt);
                            error += std::abs(f_arr(i, j, k, 1) - 0.0_rt);
                            error += std::abs(f_arr(i, j, k, 2) - 0.0_rt);
                        }
                    } else {
                        // First ghost cells only
                        const bool is_ibdy =
                            ((i == -1 || i == 8) && (j >= 0 && j <= 7) &&
                             (k >= 0 && k <= 7));
                        const bool is_jbdy =
                            ((j == -1 || j == 8) && (i >= 0 && i <= 7) &&
                             (k >= 0 && k <= 7));
                        const bool is_kbdy =
                            ((k == -1 || k == 8) && (i >= 0 && i <= 7) &&
                             (j >= 0 && j <= 7));
                        if (is_ibdy || is_jbdy || is_kbdy) {
                            error += std::abs(f_arr(i, j, k, 0) - 1.0_rt);
                            error += std::abs(f_arr(i, j, k, 1) - 2.0_rt);
                            error += std::abs(f_arr(i, j, k, 2) - 3.0_rt);
                        } else {
                            error += std::abs(f_arr(i, j, k, 0) - 0.0_rt);
                            error += std::abs(f_arr(i, j, k, 1) - 0.0_rt);
                            error += std::abs(f_arr(i, j, k, 2) - 0.0_rt);
                        }
                    }
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}
} // namespace

class FieldFillPatchTest : public MeshTest
{
public:
    void set_up_fields()
    {
        auto& frepo = mesh().field_repo();
        m_vel = &frepo.declare_field("velocity", 3, 2, 1);
        auto mac_vels = frepo.declare_face_normal_field(
            {"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
        m_umac = mac_vels[0];
        m_vmac = mac_vels[1];
        m_wmac = mac_vels[2];

        (*m_vel).setVal(0.);
        (*m_umac).setVal(0.);
        (*m_vmac).setVal(0.);
        (*m_wmac).setVal(0.);
    }

    void set_up_dirichlet()
    {
        auto& ibctype = (*m_vel).bc_type();
        auto& fbcrec = (*m_vel).bcrec();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (int i = 0; i < (*m_vel).num_comp(); ++i) {
                fbcrec[i].setLo(dir, amrex::BCType::ext_dir);
                fbcrec[i].setHi(dir, amrex::BCType::ext_dir);
            }
        }
        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            ibctype[ori] = BC::mass_inflow;
        }
        using InflowOp =
            amr_wind::BCOpCreator<TestProfile, amr_wind::ConstDirichlet>;
        AMREX_ALWAYS_ASSERT(TestProfile(*m_vel).identifier() == "TestProfile");
        (*m_vel).register_fill_patch_op<amr_wind::FieldFillPatchOps<InflowOp>>(
            mesh(), time(), InflowOp(*m_vel));
        (*m_vel).copy_bc_to_device();
        EXPECT_TRUE((*m_vel).bc_initialized());
    }

    void prep_test()
    {
        // Default dimensions are n_cell = 8 x 8 x 8
        MeshTest::populate_parameters();
        {
            amrex::Vector<int> periodic{{0, 0, 0}};
            amrex::ParmParse pp("geometry");
            pp.addarr("is_periodic", periodic);
        }
        initialize_mesh();
        set_up_fields();
        set_up_dirichlet();
    }

    amr_wind::Field* m_vel;
    amr_wind::Field* m_umac;
    amr_wind::Field* m_vmac;
    amr_wind::Field* m_wmac;
};

TEST_F(FieldFillPatchTest, dirichlet_fp)
{
    prep_test();

    // Test fillpatch and check ghost cells
    (*m_vel).fillpatch(sim().time().current_time());
    const auto err = get_field_err(*m_vel, true);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(FieldFillPatchTest, dirichlet_inflow)
{
    prep_test();

    // Test set_inflow and check boundary cells
    amr_wind::nodal_projection::set_inflow_velocity(
        sim().physics_manager(), (*m_vel), 0, time().current_time(),
        (*m_vel)(0), 1);
    const auto err = get_field_err(*m_vel, false);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(FieldFillPatchTest, dirichlet_sibling_fp)
{
    prep_test();

    // Test sibling fillpatch and check ghost cells
    amrex::Array<amr_wind::Field*, AMREX_SPACEDIM> mac_vel = {
        AMREX_D_DECL(m_umac, m_vmac, m_wmac)};
    (*m_vel).fillpatch_sibling_fields(
        time().current_time(), (*m_umac).num_grow(), mac_vel);
    auto err = get_field_err(*m_umac, true, 0);
    EXPECT_DOUBLE_EQ(err, 0.);
    err = get_field_err(*m_vmac, true, 1);
    EXPECT_DOUBLE_EQ(err, 0.);
    err = get_field_err(*m_wmac, true, 2);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(FieldFillPatchTest, dirichlet_sibling_inflow)
{
    prep_test();

    // Test sibling set_inflow and check ghost cells
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mac_vel = {
        AMREX_D_DECL(&(*m_umac)(0), &(*m_vmac)(0), &(*m_wmac)(0))};
    (*m_vel).set_inflow_sibling_fields(0, time().current_time(), mac_vel);
    auto err = get_field_err(*m_umac, false, 0);
    EXPECT_DOUBLE_EQ(err, 0.);
    err = get_field_err(*m_vmac, false, 1);
    EXPECT_DOUBLE_EQ(err, 0.);
    err = get_field_err(*m_wmac, false, 2);
    EXPECT_DOUBLE_EQ(err, 0.);
}

} // namespace amr_wind_tests
