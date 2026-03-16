#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

#include "amr-wind/utilities/FieldPlaneAveragingFine.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class FieldPlaneAveragingFineTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{8.0_rt, 8.0_rt, 8.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("amr");
            const amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 1);
            pp.add("max_grid_size", m_nx);
            pp.add("blocking_factor", 2);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{0, 0, 0}};
            pp.addarr("is_periodic", periodic);
        }
        // Create the "input file"
        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "0 0 " << z_fine_lo_in << " 8 8 " << z_fine_hi_in << '\n';

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
            new amr_wind::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    const int m_nx{8};

public:
    void test_dir(int /*dir*/);
    const amrex::Real z_fine_lo = 2;
    const amrex::Real z_fine_hi = 4;
    const amrex::Real z_fine_lo_in = z_fine_lo + 0.1_rt;
    const amrex::Real z_fine_hi_in = z_fine_hi - 0.1_rt;
};

namespace {
void init_field_linear(
    amr_wind::Field& fld,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> a,
    int dir)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real off = 0.5_rt;
                // Set values at boundaries if outside domain
                if (k < 0) {
                    off += 0.5_rt;
                }
                if (k > (8 * (lev + 1)) - 1) {
                    off -= 0.5_rt;
                }
                const amrex::GpuArray<amrex::Real, 3> x = {
                    problo[0] + ((i + 0.5_rt) * dx[0]),
                    problo[1] + ((j + 0.5_rt) * dx[1]),
                    problo[2] + ((k + off) * dx[2])};
                farrs[nbx](i, j, k, 0) = x[dir] * a[0];
                farrs[nbx](i, j, k, 1) = x[dir] * a[1];
                farrs[nbx](i, j, k, 2) = x[dir] * a[2];
            });
    }
    amrex::Gpu::streamSynchronize();
}
} // namespace

TEST_F(FieldPlaneAveragingFineTest, test_linear_fine_only)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u0 = {
        {1.0_rt, 3.5_rt, 5.6_rt}};

    populate_parameters();
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 1);

    constexpr int dir = 2;
    init_field_linear(velocityf, u0, dir);

    amr_wind::FieldPlaneAveragingFine pa_fine(
        velocityf, sim().time(), dir, true);
    pa_fine();

    constexpr int n = 20;
    const amrex::Real L = z_fine_hi - z_fine_lo;
    const amrex::Real dz = L / ((amrex::Real)n);

    // test along a line at n equidistant points in the fine zone
    for (int i = 0; i < n; ++i) {

        const amrex::Real z = z_fine_lo + (i * dz);

        const amrex::Array<amrex::Real, 3> u = {
            pa_fine.line_average_interpolated(z, 0),
            pa_fine.line_average_interpolated(z, 1),
            pa_fine.line_average_interpolated(z, 2)};

        // test each velocity field u = u0 + u0*x
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(u0[j] * (z), u[j], tol);
        }

        const amrex::Array<amrex::Real, 3> dudx = {
            pa_fine.line_derivative_interpolated(z, 0),
            pa_fine.line_derivative_interpolated(z, 1),
            pa_fine.line_derivative_interpolated(z, 2)};

        // test each velocity field du/dx = u0
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(u0[j], dudx[j], tol);
        }
    }
}

TEST_F(FieldPlaneAveragingFineTest, test_linear)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u0 = {
        {1.0_rt, 3.5_rt, 5.6_rt}};

    populate_parameters();
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 1);

    constexpr int dir = 2;
    init_field_linear(velocityf, u0, dir);

    amr_wind::FieldPlaneAveragingFine pa_fine(
        velocityf, sim().time(), dir, true);
    pa_fine();

    constexpr int n = 20;
    const amrex::Real L = z_fine_hi - z_fine_lo;
    const amrex::Real dz = L / (static_cast<amrex::Real>(n));
    const amrex::Real half_dz_pa = 0.25_rt;
    const int n_more = static_cast<int>((8.0_rt - 2.0_rt * half_dz_pa) / dz);

    // test along a line spanning domain, from first point to last
    for (int i = 0; i < n_more; ++i) {

        const amrex::Real z = half_dz_pa + (i * dz);

        const amrex::Array<amrex::Real, 3> u = {
            pa_fine.line_average_interpolated(z, 0),
            pa_fine.line_average_interpolated(z, 1),
            pa_fine.line_average_interpolated(z, 2)};

        // test each velocity field u = u0 + u0*x
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(u0[j] * (z), u[j], tol);
        }

        const amrex::Array<amrex::Real, 3> dudx = {
            pa_fine.line_derivative_interpolated(z, 0),
            pa_fine.line_derivative_interpolated(z, 1),
            pa_fine.line_derivative_interpolated(z, 2)};

        // test each velocity field du/dx = u0
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(u0[j], dudx[j], tol);
        }
    }
}

} // namespace amr_wind_tests
