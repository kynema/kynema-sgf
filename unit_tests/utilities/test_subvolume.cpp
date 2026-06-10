#include <filesystem>
#include <limits>

#include "ks_test_utils/MeshTest.H"
#include "src/utilities/subvolume/Subvolume.H"
#include "AMReX_MFIter.H"
#include "AMReX_PlotFileUtil.H"
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
            amrex::Vector<amrex::Real> probhi{{128.0_rt, 128.0_rt, 256.0_rt}};

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
        pp.addarr("labels", amrex::Vector<std::string>{"chunk1", "chunk2"});
        pp.addarr("fields", amrex::Vector<std::string>{"density"});
    }
    {
        amrex::ParmParse pp("subvolume.chunk1");
        pp.add("type", std::string("Rectangular"));
        pp.addarr("origin", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
        pp.addarr("num_points", amrex::Vector<int>{32, 32, 64});
        pp.add("dx", 4.0_rt);
    }
    {
        amrex::ParmParse pp("subvolume.chunk2");
        pp.add("type", std::string("Rectangular"));
        pp.addarr("origin", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
        pp.addarr("num_points", amrex::Vector<int>{4, 4, 4});
        pp.addarr(
            "dx_vec", amrex::Vector<amrex::Real>{8.0_rt, 16.0_rt, 12.0_rt});
        pp.add("verbose", 1);
    }

    kynema_sgf::subvolume::Subvolume subvol(sim(), "subvolume");
    subvol.initialize();
    subvol.output_actions();

    {

        const std::filesystem::path header_path(
            "post_processing/subvolume_chunk100000/Header");
        EXPECT_TRUE(std::filesystem::exists(header_path));

        amrex::PlotFileData pf("post_processing/subvolume_chunk100000");
        EXPECT_EQ(pf.finestLevel(), 0);
        EXPECT_EQ(pf.nComp(), 1);
        ASSERT_EQ(pf.varNames().size(), 1);
        EXPECT_EQ(pf.varNames()[0], "density");

        auto mf_out = pf.get(0, "density");
        for (amrex::MFIter mfi(mf_out); mfi.isValid(); ++mfi) {
            const auto bx = mfi.validbox();
            const auto arr = mf_out.const_array(mfi);
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        const amrex::Real dx_sv = 4.0_rt;
                        const amrex::Real dx_dom = 128.0_rt / 32.0_rt;
                        const amrex::Real expected =
                            dx_sv * static_cast<amrex::Real>(i + j + k) +
                            3.0_rt * 0.5_rt * dx_dom;
                        EXPECT_NEAR(
                            arr(i, j, k), expected,
                            std::numeric_limits<amrex::Real>::epsilon() *
                                1.0e3_rt);
                    }
                }
            }
        }
    }

    {

        const std::filesystem::path header_path(
            "post_processing/subvolume_chunk200000/Header");
        EXPECT_TRUE(std::filesystem::exists(header_path));

        amrex::PlotFileData pf("post_processing/subvolume_chunk200000");
        EXPECT_EQ(pf.finestLevel(), 0);
        EXPECT_EQ(pf.nComp(), 1);
        ASSERT_EQ(pf.varNames().size(), 1);
        EXPECT_EQ(pf.varNames()[0], "density");

        auto mf_out = pf.get(0, "density");
        for (amrex::MFIter mfi(mf_out); mfi.isValid(); ++mfi) {
            const auto bx = mfi.validbox();
            const auto arr = mf_out.const_array(mfi);
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        const amrex::Real dx_sv = 8.0_rt;
                        const amrex::Real dy_sv = 16.0_rt;
                        const amrex::Real dz_sv = 12.0_rt;
                        const amrex::Real dx_dom = 128.0_rt / 32.0_rt;
                        const amrex::Real expected =
                            dx_sv * static_cast<amrex::Real>(i) +
                            dy_sv * static_cast<amrex::Real>(j) +
                            dz_sv * static_cast<amrex::Real>(k) +
                            3.0_rt * 0.5_rt * dx_dom;
                        EXPECT_NEAR(
                            arr(i, j, k), expected,
                            std::numeric_limits<amrex::Real>::epsilon() *
                                1.0e3_rt);
                    }
                }
            }
        }
    }

    std::error_code ec;
    std::filesystem::remove_all("post_processing", ec);
}

} // namespace kynema_sgf_tests
