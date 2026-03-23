#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/sampling/FreeSurfaceSampler.H"
#include "amr-wind/utilities/tagging/FieldRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_vof(amr_wind::Field& vof_fld, amrex::Real water_level)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // Since VOF is cell centered
    amrex::Real offset = 0.5_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = vof_fld(lev).arrays();

        amrex::ParallelFor(
            vof_fld(lev), vof_fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);
                const amrex::Real local_vof = amrex::min<amrex::Real>(
                    1.0_rt,
                    amrex::max<amrex::Real>(
                        0.0_rt, (water_level - (z - offset * dx[2])) / dx[2]));
                farrs[nbx](i, j, k) = local_vof;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_vof_multival(
    amr_wind::Field& vof_fld, amrex::Real wl0, amrex::Real wl1, amrex::Real wl2)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // This function initializes a vof field divided by a liquid-gas interface
    // at wl2, above which is another layer of liquid, introducing interfaces
    // at wl1 and wl0. From top down: wl0, wl1, wl2.

    // Since VOF is cell centered
    amrex::Real offset = 0.5_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = vof_fld(lev).arrays();

        amrex::ParallelFor(
            vof_fld(lev), vof_fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);
                amrex::Real local_vof;
                // Above wl1
                if (z - (offset * dx[2]) > wl1) {
                    local_vof = amrex::min<amrex::Real>(
                        1.0_rt,
                        amrex::max<amrex::Real>(
                            0.0_rt, (wl0 - (z - offset * dx[2])) / dx[2]));
                } else {
                    // Above wl0
                    if (z - (offset * dx[2]) > wl2) {
                        local_vof = amrex::min<amrex::Real>(
                            1.0_rt,
                            amrex::max<amrex::Real>(
                                0.0_rt, ((z + offset * dx[2]) - wl1) / dx[2]));
                    } else {
                        // Bottom portion
                        local_vof = amrex::min<amrex::Real>(
                            1.0_rt,
                            amrex::max<amrex::Real>(
                                0.0_rt, (wl2 - (z - offset * dx[2])) / dx[2]));
                    }
                }
                farrs[nbx](i, j, k) = local_vof;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_vof_slope(
    amr_wind::Field& vof_fld,
    amrex::Real water_level,
    amrex::Real slope,
    amrex::Real domain_length)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // The interface has a slope in x and y and intersects the center of the
    // domain at the designated water level Since VOF is cell centered
    amrex::Real offset = 0.5_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = vof_fld(lev).arrays();

        amrex::ParallelFor(
            vof_fld(lev), vof_fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + offset) * dx[0]);
                const amrex::Real y = problo[1] + ((j + offset) * dx[1]);
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                // Find height of interface at current x, y
                const amrex::Real local_ht =
                    water_level + (slope * (x - 0.5_rt * domain_length)) +
                    (slope * (y - 0.5_rt * domain_length));

                const amrex::Real local_vof = amrex::min<amrex::Real>(
                    1.0_rt,
                    amrex::max<amrex::Real>(
                        0.0_rt, (local_ht - (z - offset * dx[2])) / dx[2]));
                farrs[nbx](i, j, k) = local_vof;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_vof_diffuse(amr_wind::Field& vof_fld, amrex::Real water_level)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // Since VOF is cell centered
    amrex::Real offset = 0.5_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = vof_fld(lev).arrays();

        amrex::ParallelFor(
            vof_fld(lev), vof_fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);
                const amrex::Real local_vof = amrex::min<amrex::Real>(
                    1.0_rt, amrex::max<amrex::Real>(
                                0.0_rt, 0.5_rt + (0.25_rt * (water_level - z) /
                                                  dx[2])));
                farrs[nbx](i, j, k) = local_vof;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_vof_outliers(
    amr_wind::Field& vof_fld, amrex::Real water_level, bool distort_above)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // Since VOF is cell centered
    amrex::Real offset = 0.5_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = vof_fld(lev).arrays();

        amrex::ParallelFor(
            vof_fld(lev), vof_fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);
                if (std::abs(water_level - z) < 0.5_rt * dx[2]) {
                    farrs[nbx](i, j, k) =
                        (water_level - (z - 0.5_rt * dx[2])) / dx[2];
                } else if (z > water_level) {
                    farrs[nbx](i, j, k) = distort_above ? 0.1_rt : 0.0_rt;
                } else {
                    farrs[nbx](i, j, k) = !distort_above ? 0.9_rt : 1.0_rt;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Custom mesh class to be able to refine like a simulation would
//  - combination of AmrTestMesh and incflo classes
//  - with ability to initialize the refiner and regrid
class FSRefineMesh : public AmrTestMesh
{
public:
    FSRefineMesh() : m_mesh_refiner(new amr_wind::RefineCriteriaManager(m_sim))
    {}
    void init_refiner() { m_mesh_refiner->initialize(); }
    void remesh() { regrid(0, 0.0_rt); }

protected:
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override
    {
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        m_repo.make_new_level_from_scratch(lev, time, ba, dm);
    }

    void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override
    {
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        m_repo.make_new_level_from_coarse(lev, time, ba, dm);
    }

    void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override
    {
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        m_repo.remake_level(lev, time, ba, dm);
    }

    void ClearLevel(int lev) override { m_repo.clear_level(lev); }

    void ErrorEst(
        int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override
    {
        m_mesh_refiner->tag_cells(lev, tags, time, ngrow);
    }

private:
    std::unique_ptr<amr_wind::RefineCriteriaManager> m_mesh_refiner;
};

class FreeSurfaceImpl : public amr_wind::sampling::FreeSurfaceSampler
{
public:
    FreeSurfaceImpl(amr_wind::CFDSim& sim)
        : amr_wind::sampling::FreeSurfaceSampler(sim)
    {}
    int check_output(const std::string& op, amrex::Real check_val);
    int check_output(int cidx, const std::string& op, amrex::Real check_val);
    int check_output_vec(
        const std::string& op, amrex::Vector<amrex::Real> check_val);
    int check_pos(int cidx, const std::string& op, amrex::Real check_val);
    int check_sloc(const std::string& op);

protected:
    const amrex::Real m_tol = std::numeric_limits<float>::epsilon() * 1.0e1_rt;
};

int FreeSurfaceImpl::check_output(
    int cidx, const std::string& op, amrex::Real check_val)
{
    // Get number of points and output array
    auto npts_tot = num_gridpoints();
    auto out = heights();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(out[(cidx * npts_tot) + n], check_val);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(out[(cidx * npts_tot) + n], check_val);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(out[(cidx * npts_tot) + n], check_val, m_tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

int FreeSurfaceImpl::check_output(const std::string& op, amrex::Real check_val)
{
    return check_output(0, op, check_val);
}

int FreeSurfaceImpl::check_output_vec(
    const std::string& op, amrex::Vector<amrex::Real> check_val)
{
    // Get number of points and output array
    auto npts_tot = num_gridpoints();
    auto out = heights();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(out[n], check_val[n]);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(out[n], check_val[n]);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(out[n], check_val[n], m_tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

int FreeSurfaceImpl::check_pos(
    const int cidx, const std::string& op, amrex::Real check_val)
{
    // Get number of points and position array
    auto npts_tot = num_gridpoints();
    auto locs = grid_locations();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(locs[n][cidx], check_val);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(locs[n][cidx], check_val);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(locs[n][cidx], check_val, m_tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

int FreeSurfaceImpl::check_sloc(const std::string& op)
{
    // Get number of points and sampling locations array
    auto npts_tot = num_points();
    amr_wind::sampling::SampleLocType sample_locs;
    output_locations(sample_locs);
    // Get locations from other functions
    auto gridlocs = grid_locations();
    auto out = heights();
    // Loop through grid points and check output
    int icheck = 0;
    const auto& locs = sample_locs.locations();
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(locs[n][0], gridlocs[n][0]);
            EXPECT_EQ(locs[n][1], gridlocs[n][1]);
            EXPECT_EQ(locs[n][2], out[n]);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(locs[n][0], gridlocs[n][0]);
                EXPECT_LT(locs[n][1], gridlocs[n][1]);
                EXPECT_LT(locs[n][2], out[n]);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(locs[n][0], gridlocs[n][0], m_tol);
                    EXPECT_NEAR(locs[n][1], gridlocs[n][1], m_tol);
                    EXPECT_NEAR(locs[n][2], out[n], m_tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

} // namespace

class FreeSurfaceTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 64}};
            pp.add("max_level", m_nlev);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");

            pp.addarr("prob_lo", m_problo);
            pp.addarr("prob_hi", m_probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    void setup_grid_0d(int ninst, const std::string& fsname)
    {
        amrex::ParmParse pp(fsname);
        pp.add("output_interval", 1);
        pp.add("num_instances", ninst);
        pp.addarr("plane_num_points", amrex::Vector<int>{1, 1});
        pp.addarr("plane_start", m_pt_coord);
        pp.addarr("plane_end", m_pt_coord);
    }
    void setup_grid_0d(int ninst) { setup_grid_0d(ninst, "freesurface"); }
    void setup_grid_0d(int ninst, amrex::Real xy_sample)
    {
        amrex::ParmParse pp("freesurface");
        pp.add("output_interval", 1);
        pp.add("num_instances", ninst);
        pp.addarr("plane_num_points", amrex::Vector<int>{1, 1});
        const amrex::Vector<amrex::Real> pt_coord{
            xy_sample, xy_sample, m_pt_coord[2]};
        pp.addarr("plane_start", pt_coord);
        pp.addarr("plane_end", pt_coord);
    }
    void setup_grid_2d(int ninst)
    {
        amrex::ParmParse pp("freesurface");
        pp.add("output_interval", 1);
        pp.add("num_instances", ninst);
        pp.addarr("plane_num_points", amrex::Vector<int>{npts, npts});
        pp.addarr("plane_start", m_pl_start);
        pp.addarr("plane_end", m_pl_end);
    }
    void setup_grid_2d_narrow(const std::string& fsname)
    {
        amrex::ParmParse pp(fsname);
        pp.add("output_interval", 1);
        pp.add("num_instances", 1);
        pp.addarr("plane_num_points", amrex::Vector<int>{npts, npts});
        pp.addarr("plane_start", m_plnarrow_s);
        pp.addarr("plane_end", m_plnarrow_e);
    }
    void setup_grid_2d_narrow() { setup_grid_2d_narrow("freesurface"); }
    void setup_fieldrefinement()
    {
        amrex::ParmParse pp("tagging");
        pp.add("labels", (std::string) "t1");
        amrex::ParmParse ppt1("tagging.t1");
        ppt1.add("type", (std::string) "FieldRefinement");
        ppt1.add("field_name", m_fname);
        ppt1.addarr("field_error", amrex::Vector<amrex::Real>{m_fref_val});
    }
    // Parameters to reuse
    const amrex::Real m_water_level0 = 64.0_rt;
    const amrex::Real m_water_level1 = 31.5_rt;
    const amrex::Real m_water_level2 = 65.0_rt;
    const amrex::Vector<amrex::Real> m_problo{{0.0_rt, 0.0_rt, -4.0_rt}};
    const amrex::Vector<amrex::Real> m_probhi{{128.0_rt, 128.0_rt, 124.0_rt}};
    const amrex::Vector<amrex::Real> m_pt_coord{{63.0_rt, 65.0_rt, 0.0_rt}};
    const amrex::Vector<amrex::Real> m_pl_start{{0.0_rt, 0.0_rt, 0.0_rt}};
    const amrex::Vector<amrex::Real> m_pl_end{{128.0_rt, 128.0_rt, 0.0_rt}};
    const amrex::Vector<amrex::Real> m_plnarrow_s{{63.0_rt, 63.0_rt, 0.0_rt}};
    const amrex::Vector<amrex::Real> m_plnarrow_e{{65.0_rt, 65.0_rt, 0.0_rt}};
    static constexpr int npts = 3;
    const amrex::Real m_fref_val = 0.5_rt;
    const std::string m_fname = "flag";
    int m_nlev = 0;
};

TEST_F(FreeSurfaceTest, point)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid_0d(1);

    init_vof(vof, m_water_level0);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Check number of points
    auto ngp = tool.num_gridpoints();
    EXPECT_EQ(ngp, 1);
    // Check location after being read
    int npos = tool.check_pos(0, "=", m_pt_coord[0]);
    ASSERT_EQ(npos, 1);
    npos = tool.check_pos(1, "=", m_pt_coord[1]);
    ASSERT_EQ(npos, 1);
    // Check output value
    int nout = tool.check_output("~", m_water_level0);
    ASSERT_EQ(nout, 1);
    // Check sampling locations
    int nsloc = tool.check_sloc("~");
    ASSERT_EQ(nsloc, 1);
}

TEST_F(FreeSurfaceTest, plane)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid_2d(1);

    init_vof(vof, m_water_level1);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Check number of points
    auto ngp = tool.num_gridpoints();
    EXPECT_EQ(ngp, npts * npts);
    // Check output value
    int nout = tool.check_output("~", m_water_level1);
    ASSERT_EQ(nout, npts * npts);
}

TEST_F(FreeSurfaceTest, multivalued)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid_2d(3);

    // Parameters of water level
    amrex::Real wl0 = m_probhi[2] * 2 / 3;
    amrex::Real wl1 = m_probhi[2] / 2;
    amrex::Real wl2 = m_probhi[2] / 5;

    init_vof_multival(vof, wl0, wl1, wl2);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Check number of outputs
    auto heights = tool.heights();
    EXPECT_EQ(heights.size(), 3 * npts * npts);

    // Check output values
    int nout = tool.check_output(0, "~", wl0);
    ASSERT_EQ(nout, npts * npts);
    nout = tool.check_output(1, "~", wl1);
    ASSERT_EQ(nout, npts * npts);
    nout = tool.check_output(2, "~", wl2);
    ASSERT_EQ(nout, npts * npts);
}

TEST_F(FreeSurfaceTest, sloped)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid_2d_narrow();

    amrex::Real slope = 0.125_rt;
    amrex::Real domain_l = m_probhi[0];
    init_vof_slope(vof, m_water_level2, slope, domain_l);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Calculate expected output values
    amrex::Vector<amrex::Real> out_vec(static_cast<long>(npts * npts), 0.0_rt);
    // Step in x, then y
    out_vec[0] = (m_water_level2 + (slope * (-1.0_rt - 1.0_rt)));
    out_vec[1] = (m_water_level2 + (slope * (+0.0_rt - 1.0_rt)));
    out_vec[2] = (m_water_level2 + (slope * (+1.0_rt - 1.0_rt)));
    out_vec[3] = (m_water_level2 + (slope * (-1.0_rt + 0.0_rt)));
    out_vec[4] = (m_water_level2 + (slope * (+0.0_rt + 0.0_rt)));
    out_vec[5] = (m_water_level2 + (slope * (+1.0_rt + 0.0_rt)));
    out_vec[6] = (m_water_level2 + (slope * (-1.0_rt + 1.0_rt)));
    out_vec[7] = (m_water_level2 + (slope * (+0.0_rt + 1.0_rt)));
    out_vec[8] = (m_water_level2 + (slope * (+1.0_rt + 1.0_rt)));
    // Check output value
    int nout = tool.check_output_vec("~", out_vec);
    ASSERT_EQ(nout, npts * npts);
}

TEST_F(FreeSurfaceTest, multisampler)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);

    // Set up parameters for one sampler
    setup_grid_0d(1, "freesurface0");

    // Set up parameters for another sampler
    setup_grid_2d_narrow("freesurface1");

    // Initialize VOF distribution and access sim
    init_vof(vof, m_water_level1);
    auto& m_sim = sim();

    // Initialize first sampler
    FreeSurfaceImpl tool1(m_sim);
    // Populate label, would be done by Sampling
    tool1.label() = "0";
    tool1.initialize("freesurface0");

    // Initialize second sampler
    FreeSurfaceImpl tool2(m_sim);
    tool1.label() = "1";
    tool2.initialize("freesurface1");
}

TEST_F(FreeSurfaceTest, regrid)
{
    // Allow 2 levels
    m_nlev = 1;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_grid_2d(1);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FSRefineMesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0_rt * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize VOF distribution and access sim
    init_vof(vof, m_water_level1);
    auto& rsim = rmesh.sim();

    // Initialize sampler and check result on initial mesh
    FreeSurfaceImpl tool(rsim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();
    tool.check_output("~", m_water_level1);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0);

    // Regrid, update fields, and do post_regrid_actions for tool
    rmesh.remesh();
    tool.post_regrid_actions();

    // Check that result is unchanged on new mesh
    tool.update_sampling_locations();
    tool.check_output("~", m_water_level1);
}

TEST_F(FreeSurfaceTest, point_diffuse)
{
    initialize_mesh();
    sim().activate_overset();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    auto& iblank = repo.get_int_field("iblank_cell");
    iblank.setVal(-1);
    setup_grid_0d(1);

    // -- Chosen to be cell center --
    amrex::Real water_lev_diffuse = 61.;
    init_vof_diffuse(vof, water_lev_diffuse);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Check number of points
    auto ngp = tool.num_gridpoints();
    EXPECT_EQ(ngp, 1);
    // Check location after being read
    int npos = tool.check_pos(0, "=", m_pt_coord[0]);
    ASSERT_EQ(npos, 1);
    npos = tool.check_pos(1, "=", m_pt_coord[1]);
    ASSERT_EQ(npos, 1);
    // Check output value
    int nout = tool.check_output("~", water_lev_diffuse);
    ASSERT_EQ(nout, 1);
    // Check sampling locations
    int nsloc = tool.check_sloc("~");
    ASSERT_EQ(nsloc, 1);

    // -- Chosen to be cell edge --
    water_lev_diffuse = 62.;
    init_vof_diffuse(vof, water_lev_diffuse);
    tool.update_sampling_locations();

    // Check output value
    nout = tool.check_output("~", water_lev_diffuse);
    ASSERT_EQ(nout, 1);

    // -- Chosen to be neither cell center nor edge --
    water_lev_diffuse = 61.5_rt;
    init_vof_diffuse(vof, water_lev_diffuse);
    tool.update_sampling_locations();

    // Check output value
    nout = tool.check_output("~", water_lev_diffuse);
    ASSERT_EQ(nout, 1);
}

TEST_F(FreeSurfaceTest, point_outliers)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid_0d(1);
    {
        amrex::ParmParse pp("freesurface");
        pp.add("linear_interp_extent_from_xhi", 66.1_rt);
    }

    amrex::Real water_lev = 61.5_rt;
    // Sets up field with multiphase cells above interface
    init_vof_outliers(vof, water_lev, true);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Linear interpolation will get different value
    amrex::Real local_vof = (water_lev - 60.0_rt) / 2.0_rt;
    amrex::Real water_lev_linear =
        61.0_rt + ((0.5_rt - local_vof) * (2.0_rt / (0.1_rt - local_vof)));

    // Check output value
    int nout = tool.check_output("~", water_lev_linear);
    ASSERT_EQ(nout, 1);
    // Check sampling locations
    int nsloc = tool.check_sloc("~");
    ASSERT_EQ(nsloc, 1);

    // Now distort VOF field below interface
    init_vof_outliers(vof, water_lev, false);
    tool.update_sampling_locations();
    water_lev_linear =
        61.0_rt + ((0.5_rt - local_vof) * (2.0_rt / (0.0_rt - local_vof)));
    nout = tool.check_output("~", water_lev_linear);
    ASSERT_EQ(nout, 1);
    nsloc = tool.check_sloc("~");
    ASSERT_EQ(nsloc, 1);

    // Repeat with target cell having 0 < vof < 0.5_rt
    water_lev = 60.5_rt;
    init_vof_outliers(vof, water_lev, true);
    tool.update_sampling_locations();
    local_vof = (water_lev - 60.0_rt) / 2.0_rt;
    water_lev_linear =
        61.0_rt - ((0.5_rt - local_vof) * (2.0_rt / (1.0_rt - local_vof)));
    nout = tool.check_output("~", water_lev_linear);
    ASSERT_EQ(nout, 1);
    init_vof_outliers(vof, water_lev, false);
    tool.update_sampling_locations();
    water_lev_linear =
        61.0_rt - ((0.5_rt - local_vof) * (2.0_rt / (0.9_rt - local_vof)));
    nout = tool.check_output("~", water_lev_linear);
    ASSERT_EQ(nout, 1);
    nsloc = tool.check_sloc("~");
    ASSERT_EQ(nsloc, 1);
}

TEST_F(FreeSurfaceTest, point_outliers_off)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid_0d(1);
    {
        amrex::ParmParse pp("freesurface");
        pp.add("linear_interp_extent_from_xhi", 64);
    }

    amrex::Real water_lev = 61.5_rt;
    // Sets up field with multiphase cells above interface
    init_vof_outliers(vof, water_lev, true);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Because linear interpolation is off and there are multiphase cells all
    // the way above the interface, it will get a bad value, middle of top cell
    const amrex::Real water_lev_geom = 123.;

    // Check output value
    int nout = tool.check_output("~", water_lev_geom);
    ASSERT_EQ(nout, 1);
}

TEST_F(FreeSurfaceTest, point_diffuse_in_single_phase)
{
    initialize_mesh();
    sim().activate_overset();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    auto& iblank = repo.get_int_field("iblank_cell");
    iblank.setVal(-1);
    // Set up sampler to be on lateral edge of cell
    setup_grid_0d(1, 64.0_rt - 0.1_rt);

    // Set up water level to be just below bottom edge of cell
    const amrex::Real water_lev_diffuse = 60.0_rt - 0.1_rt;
    const amrex::Real vof_slope = 0.1_rt;
    // Set up vof with nonzero slope but 0 in current cell
    init_vof_slope(vof, water_lev_diffuse, vof_slope, 124.0_rt);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim);
    tool.initialize("freesurface");
    tool.update_sampling_locations();

    // Check output value
    const amrex::Real vof_cell = 0.0_rt;
    const amrex::Real vof_mz = (water_lev_diffuse - 58.0_rt) / 2.0_rt;
    const amrex::Real vof_pxy =
        (water_lev_diffuse + 4.0_rt * vof_slope - 60.0_rt) / 2.0_rt;
    const amrex::Real vof_c =
        vof_cell + (2.0_rt * (vof_pxy - vof_cell) / 4.0_rt * (2.0_rt - 0.1_rt));
    const amrex::Real slope_z = (vof_cell - vof_mz) / 2.0_rt;
    const amrex::Real ht =
        61.0_rt + ((0.5_rt - vof_c) / (slope_z + amr_wind::constants::EPS));
    const amrex::Real ht_est =
        water_lev_diffuse + (2.0_rt * vof_slope * (2.0_rt - 0.1_rt));
    int nout = tool.check_output("~", ht);
    ASSERT_EQ(nout, 1);
    // The linear interpolation answer should be close to the naive calc
    EXPECT_NEAR(ht, ht_est, 5.0e-2_rt);
}

} // namespace amr_wind_tests
