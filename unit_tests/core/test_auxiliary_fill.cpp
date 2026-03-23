#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void auxiliary_fill_boundary(
    amr_wind::Field& velocity, amr_wind::IntField& indices, const int comp = 0)
{
    const int nlevels = velocity.repo().num_active_levels();
    const int ncomp = velocity.num_comp();

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const amrex::Box& minBox =
                velocity.repo().mesh().boxArray(lev).minimalBox();
            amrex::IntVect plo(minBox.loVect());
            amrex::IntVect phi(minBox.hiVect());
            const int normal = ori.coordDir();
            plo[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            const amrex::Box domain_bdy_bx(plo, phi);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(velocity(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {

                auto sbx = mfi.growntilebox(1);
                auto shift_to_cc = amrex::IntVect(0);
                const auto& bx =
                    amr_wind::utils::face_aware_boundary_box_intersection(
                        shift_to_cc, sbx, domain_bdy_bx, ori);
                if (bx.isEmpty()) {
                    continue;
                }

                const auto& dest = velocity(lev).array(mfi);
                const auto& idx = indices(lev).array(mfi);
                amrex::ParallelFor(
                    bx, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        dest(i, j, k, n) =
                            static_cast<amrex::Real>(comp + n + 1);
                        idx(i + shift_to_cc[0], j + shift_to_cc[1],
                            k + shift_to_cc[2], (3 * comp) + (3 * n)) =
                            shift_to_cc[0];
                        idx(i + shift_to_cc[0], j + shift_to_cc[1],
                            k + shift_to_cc[2], (3 * comp) + (3 * n) + 1) =
                            shift_to_cc[1];
                        idx(i + shift_to_cc[0], j + shift_to_cc[1],
                            k + shift_to_cc[2], (3 * comp) + (3 * n) + 2) =
                            shift_to_cc[2];
                    });
            }
        }
    }
}

amrex::Real get_field_err(
    amr_wind::Field& field,
    amr_wind::IntField& indices,
    const bool check_all_ghost,
    const int comp = 0)
{
    const int lev = 0;
    amrex::Real error_total = 0;
    const int ncomp = field.num_comp();

    error_total += amrex::ReduceSum(
        field(lev), indices(lev), field(lev).nGrow(),
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& f_arr,
            amrex::Array4<int const> const& i_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) {
                // Do indices manually to check against box operations in code
                if (ncomp == 1) {
                    // Checking a MAC velocity
                    const bool in_x_bdy = (i == -1 || i == 8);
                    const bool not_in_yz_bdy =
                        (j >= 0 && j <= 7) && (k >= 0 && k <= 7);
                    const bool in_y_bdy = (j == -1 || j == 8);
                    const bool not_in_xz_bdy =
                        (i >= 0 && i <= 7) && (k >= 0 && k <= 7);
                    const bool in_z_bdy = (k == -1 || k == 8);
                    const bool not_in_xy_bdy =
                        (j >= 0 && j <= 7) && (i >= 0 && i <= 7);
                    if (comp == 0) {
                        // Checking u
                        if (check_all_ghost) {
                            // Ghost cells in boundaries, but not corners
                            const bool in_xf_bdy = (i == -1 || i == 9);
                            if ((in_xf_bdy && not_in_yz_bdy) ||
                                (in_y_bdy && not_in_xz_bdy) ||
                                (in_z_bdy && not_in_xy_bdy)) {
                                error += std::abs(f_arr(i, j, k) - 1.0_rt);
                                if (i == 9) {
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i - 1, j, k, 0) +
                                                      1));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i - 1, j, k, 1) -
                                                      0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i - 1, j, k, 2) -
                                                      0));
                                } else if (i == -1) {
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 0) - 0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 1) - 0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 2) - 0));
                                }
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
                            // Ghost cells in boundaries, but not corners
                            const bool in_yf_bdy = (j == -1 || j == 9);
                            if ((in_x_bdy && not_in_yz_bdy) ||
                                (in_yf_bdy && not_in_xz_bdy) ||
                                (in_z_bdy && not_in_xy_bdy)) {
                                error += std::abs(f_arr(i, j, k) - 2.0_rt);
                                if (j == 9) {
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j - 1, k, 3) -
                                                      0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j - 1, k, 4) +
                                                      1));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j - 1, k, 5) -
                                                      0));
                                } else if (j == -1) {
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 3) - 0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 4) - 0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 5) - 0));
                                }
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
                            // Ghost cells in boundaries, but not corners
                            const bool in_zf_bdy = (k == -1 || k == 9);
                            if ((in_x_bdy && not_in_yz_bdy) ||
                                (in_y_bdy && not_in_xz_bdy) ||
                                (in_zf_bdy && not_in_xy_bdy)) {
                                error += std::abs(f_arr(i, j, k) - 3.0_rt);
                                if (k == 9) {
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k - 1, 6) -
                                                      0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k - 1, 7) -
                                                      0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k - 1, 8) +
                                                      1));
                                } else if (k == -1) {
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 6) - 0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 7) - 0));
                                    error += std::abs(
                                        (amrex::Real)(i_arr(i, j, k, 8) - 0));
                                }
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
                            for (int ni = 0; ni < 9; ++ni) {
                                error += std::abs(
                                    (amrex::Real)(i_arr(i, j, k, ni) - 0));
                            }
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

class AuxiliaryFillTest : public MeshTest
{
public:
    void set_up_fields_cc()
    {
        auto& frepo = mesh().field_repo();
        m_vel = &frepo.declare_field("velocity", 3, 1, 1);
        m_ind = &frepo.declare_int_field("indices", 9, 1);

        (*m_vel).setVal(0.);
        (*m_ind).setVal(5);
    }

    void set_up_fields_face()
    {
        auto& frepo = mesh().field_repo();
        auto mac_vels = frepo.declare_face_normal_field(
            {"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
        m_umac = mac_vels[0];
        m_vmac = mac_vels[1];
        m_wmac = mac_vels[2];
        m_ind = &frepo.declare_int_field("indices", 9, 1);

        (*m_umac).setVal(0.);
        (*m_vmac).setVal(0.);
        (*m_wmac).setVal(0.);
        (*m_ind).setVal(5);
    }

    void prep_test(const bool is_cc)
    {
        // Default dimensions are n_cell = 8 x 8 x 8
        MeshTest::populate_parameters();
        {
            amrex::Vector<int> periodic{{0, 0, 0}};
            amrex::ParmParse pp("geometry");
            pp.addarr("is_periodic", periodic);
        }
        initialize_mesh();
        if (is_cc) {
            set_up_fields_cc();
        } else {
            set_up_fields_face();
        }
    }

    amr_wind::Field* m_vel;
    amr_wind::Field* m_umac;
    amr_wind::Field* m_vmac;
    amr_wind::Field* m_wmac;
    amr_wind::IntField* m_ind;
};

TEST_F(AuxiliaryFillTest, velocity_cc)
{
    prep_test(true);

    // Do fill and check ghost cells
    auxiliary_fill_boundary(*m_vel, *m_ind);
    const auto err = get_field_err(*m_vel, *m_ind, false);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(AuxiliaryFillTest, velocity_face)
{
    prep_test(false);

    // Do fill and check ghost cells
    auxiliary_fill_boundary(*m_umac, *m_ind, 0);
    const auto u_err = get_field_err(*m_umac, *m_ind, true, 0);
    EXPECT_DOUBLE_EQ(u_err, 0.);

    auxiliary_fill_boundary(*m_vmac, *m_ind, 1);
    const auto v_err = get_field_err(*m_vmac, *m_ind, true, 1);
    EXPECT_DOUBLE_EQ(v_err, 0.);

    auxiliary_fill_boundary(*m_wmac, *m_ind, 2);
    const auto w_err = get_field_err(*m_wmac, *m_ind, true, 2);
    EXPECT_DOUBLE_EQ(w_err, 0.);
}

} // namespace amr_wind_tests
