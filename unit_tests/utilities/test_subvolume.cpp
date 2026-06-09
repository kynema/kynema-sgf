#include <filesystem>

#include "ks_test_utils/MeshTest.H"
#include "src/utilities/subvolume/Subvolume.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void init_field(kynema_sgf::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == kynema_sgf::FieldLoc::CELL) {
        offset = 0.5_rt;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(), ncomp,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                const amrex::Real x = problo[0] + ((i + offset) * dx[0]);
                const amrex::Real y = problo[1] + ((j + offset) * dx[1]);
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);
                farrs[nbx](i, j, k, n) = x + y + z;
            });
    }
    amrex::Gpu::streamSynchronize();
}

class SubvolumeTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 64}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{128.0_rt, 128.0_rt, 128.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

} // namespace

TEST_F(SubvolumeTest, rectangular_subvolume_output)
{
    initialize_mesh();

    auto& repo = sim().repo();
    auto& rho = repo.declare_field("density", 1, 2);
    init_field(rho);

    {
        amrex::ParmParse pp("subvolume");
        pp.add("output_interval", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"chunk1"});
        pp.addarr("fields", amrex::Vector<std::string>{"density"});
    }
    {
        amrex::ParmParse pp("subvolume.chunk1");
        pp.add("type", std::string("Rectangular"));
        pp.addarr(
            "origin", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
        pp.addarr("num_points", amrex::Vector<int>{4, 4, 4});
        pp.add("dx", 4.0_rt);
    }

    kynema_sgf::subvolume::Subvolume subvol(sim(), "subvolume");
    subvol.initialize();
    subvol.output_actions();

    const std::filesystem::path header_path(
        "post_processing/subvolume_chunk100000/Header");
    EXPECT_TRUE(std::filesystem::exists(header_path));

    std::error_code ec;
    std::filesystem::remove_all("post_processing", ec);
}

} // namespace kynema_sgf_tests
