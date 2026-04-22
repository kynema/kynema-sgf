#include "ks_test_utils/MeshTest.H"
#include "src/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

class FieldOpsTest : public MeshTest
{
public:
    void declare_default_fields()
    {
        auto& frepo = mesh().field_repo();
        frepo.declare_field("scalar_field", 1, 1, 2);
        frepo.declare_field("vector_field", 3, 1, 2);
    }

    static void initialise_default_fields(
        kynema_sgf::Field& field,
        amrex::Vector<amrex::Geometry> geom,
        const int nlevels)
    {
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const auto& farrs = field(lev).arrays();
            amrex::ParallelFor(
                field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                    farrs[nbx](i, j, k) = 1.0_rt - (x + y + z);
                });
        }
        amrex::Gpu::streamSynchronize();
    }
};

TEST_F(FieldOpsTest, compute_max_magnitude)
{
    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& field = frepo.declare_field("scalar_field", 3, 0, 2);
    const auto& geom = mesh().Geom();
    const int nlevels = mesh().finestLevel() + 1;

    initialise_default_fields(field, geom, nlevels);

    amrex::Real global_maximum =
        kynema_sgf::field_ops::global_max_magnitude(field);
    EXPECT_NEAR(
        global_maximum, 21.5_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

} // namespace kynema_sgf_tests
