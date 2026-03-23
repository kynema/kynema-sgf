#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "amr-wind/utilities/tagging/FieldRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_velocity(
    amr_wind::Field& vel_fld,
    amrex::Real u,
    amrex::Real v,
    amrex::Real w,
    amrex::Real lev0_fac)
{
    const int nlevels = vel_fld.repo().num_active_levels();
    const amrex::GpuArray<amrex::Real, 3> vels = {u, v, w};
    for (int lev = 0; lev < nlevels; ++lev) {

        // Apply a factor to the base level values
        const amrex::Real fac = (lev == 0) ? lev0_fac : 1.0_rt;
        const auto& farrs = vel_fld(lev).arrays();

        amrex::ParallelFor(
            vel_fld(lev), vel_fld.num_grow(), vel_fld.num_comp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                // Mix positive and negative as check on L2 norm
                farrs[nbx](i, j, k, n) =
                    fac * (i % 2 == 0 ? vels[n] : -vels[n]);
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_velocity(
    const amrex::Real var_f,
    amr_wind::Field& vel_fld,
    const amrex::Real u,
    const amrex::Real v,
    const amrex::Real w)
{
    const int nlevels = vel_fld.repo().num_active_levels();
    const amrex::GpuArray<amrex::Real, 3> vels = {u, v, w};
    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& farrs = vel_fld(lev).arrays();

        amrex::ParallelFor(
            vel_fld(lev), vel_fld.num_grow(), vel_fld.num_comp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                // Mix positive and negative as check on L2 norm
                farrs[nbx](i, j, k, n) =
                    (i % 2 == 0 ? 1.0_rt - var_f : 1.0_rt + var_f) * vels[n];
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Custom mesh class to be able to refine like a simulation would
//  - combination of AmrTestMesh and incflo classes
//  - with ability to initialize the refiner and regrid
class FNRefinemesh : public AmrTestMesh
{
public:
    FNRefinemesh() : m_mesh_refiner(new amr_wind::RefineCriteriaManager(m_sim))
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

class FieldNormsImpl : public amr_wind::field_norms::FieldNorms
{
public:
    FieldNormsImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::field_norms::FieldNorms(sim, label)
    {}
    void check_output(
        amrex::Real check_val0, amrex::Real check_val1, amrex::Real check_val2);
    void check_output(amrex::Real check_val);
    void check_output4(amrex::Real check_val);

protected:
    // No file output during test
    void prepare_ascii_file() override {}
    void write_ascii() override {}
    const amrex::Real m_tol = std::numeric_limits<float>::epsilon() * 1.0e1_rt;
};

void FieldNormsImpl::check_output(
    amrex::Real check_val0, amrex::Real check_val1, amrex::Real check_val2)
{
    // Get variable names and check
    EXPECT_EQ(var_names()[0], (std::string) "velocityx");
    EXPECT_EQ(var_names()[1], (std::string) "velocityy");
    EXPECT_EQ(var_names()[2], (std::string) "velocityz");
    // Loop through norm values and check them
    const amrex::Real tol = amr_wind::constants::TIGHT_TOL;
    EXPECT_NEAR(field_norms()[0], check_val0, tol);
    EXPECT_NEAR(field_norms()[1], check_val1, tol);
    EXPECT_NEAR(field_norms()[2], check_val2, tol);
}

void FieldNormsImpl::check_output(amrex::Real check_val)
{
    // Get variable names and check
    EXPECT_EQ(var_names()[0], (std::string) "velocity");
    EXPECT_NEAR(field_norms()[0], check_val, amr_wind::constants::TIGHT_TOL);
}

void FieldNormsImpl::check_output4(amrex::Real check_val)
{
    // Get variable names and check
    EXPECT_EQ(var_names()[0], (std::string) "pressure");
    EXPECT_EQ(var_names()[1], (std::string) "u_mac");
    EXPECT_EQ(var_names()[2], (std::string) "v_mac");
    EXPECT_EQ(var_names()[3], (std::string) "w_mac");
    // Loop through norm values and check them
    const amrex::Real tol = amr_wind::constants::TIGHT_TOL;
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(field_norms()[n], check_val, tol * check_val * 0.1_rt);
    }
}

} // namespace

class FieldNormsTest : public MeshTest
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
    static void setup_fnorm(bool levelmask_flag)
    {
        amrex::ParmParse pp("fieldnorm");
        pp.add("output_interval", 1);
        if (!levelmask_flag) {
            pp.add("mask_redundant_grids", false);
        }
    }
    void setup_fieldrefinement()
    {
        amrex::ParmParse pp("tagging");
        pp.add("labels", (std::string) "t1");
        amrex::ParmParse ppt1("tagging.t1");
        ppt1.add("type", (std::string) "FieldRefinement");
        ppt1.add("field_name", m_fname);
        ppt1.addarr("field_error", amrex::Vector<amrex::Real>{m_fref_val});
    }
    void setup_intfieldrefinement()
    {
        amrex::ParmParse pp("tagging");
        pp.add("labels", (std::string) "t1");
        amrex::ParmParse ppt1("tagging.t1");
        ppt1.add("type", (std::string) "FieldRefinement");
        ppt1.add("field_name", m_ifname);
        ppt1.addarr("field_error", amrex::Vector<int>{m_ifref_val});
    }
    static void setup_tagging_box()
    {
        amrex::ParmParse ppt1("tagging.t1");
        amrex::Vector<amrex::Real> blo = {3.0_rt, 3.0_rt, 3.0_rt};
        amrex::Vector<amrex::Real> bhi = {30.0_rt, 30.0_rt, 30.0_rt};
        ppt1.addarr("box_lo", blo);
        ppt1.addarr("box_hi", bhi);
    }
    // Parameters to reuse
    const amrex::Vector<amrex::Real> m_problo{{0.0_rt, 0.0_rt, -4.0_rt}};
    const amrex::Vector<amrex::Real> m_probhi{{128.0_rt, 128.0_rt, 124.0_rt}};
    const amrex::Real m_fref_val = 0.5_rt;
    const std::string m_fname = "flag";
    const int m_ifref_val = 1;
    const std::string m_ifname = "iflag";
    const int m_nlev = 1;
    const int m_nx = 32;
    const int m_ny = 32;
    const int m_nz = 64;
    // Number of cells on base level
    const int m_ncell0 = m_nx * m_ny * m_nz;
    // Cell volume on base level
    const amrex::Real m_cv0 =
        (128.0_rt / m_nx) * (128.0_rt / m_ny) * (128.0_rt / m_nz);
    // Total domain volume
    const amrex::Real m_dv = 128.0_rt * 128.0_rt * 128.0_rt;
    // Velocity component values
    const amrex::Real m_u = 0.9_rt;
    const amrex::Real m_v = 1.3_rt;
    const amrex::Real m_w = 2.7_rt;
};

TEST_F(FieldNormsTest, levelmask_on)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0_rt * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5_rt;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.output_actions();
    // Only highest level will be counted
    // sqrt(number of cells * cell volume * cell value * cell value)
    amrex::Real unorm =
        std::sqrt((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_u * m_u / m_dv);
    amrex::Real vnorm =
        std::sqrt((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_v * m_v / m_dv);
    amrex::Real wnorm =
        std::sqrt((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.output_actions();
    // Only base level exists, which includes factor
    unorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u / m_dv);
    vnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v / m_dv);
    wnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);
}

TEST_F(FieldNormsTest, levelmask_on_with_box)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    setup_tagging_box();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0_rt * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5_rt;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.output_actions();
    tool.check_output(
        1.342655804506534_rt, 1.9393917176204616_rt, 4.0279674135195087_rt);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.output_actions();
    tool.check_output(
        1.3500000000000343_rt, 1.9499999999999589_rt, 4.0500000000000043_rt);
}

// Not sure why this option would be used, but it is available
TEST_F(FieldNormsTest, levelmask_off)
{
    bool levelmask = false;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2 * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5_rt;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.output_actions();
    // Both levels will be counted
    amrex::Real unorm = std::sqrt(
        ((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_u * m_u +
         m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u) /
        m_dv);
    amrex::Real vnorm = std::sqrt(
        ((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_v * m_v +
         m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v) /
        m_dv);
    amrex::Real wnorm = std::sqrt(
        ((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_w * m_w +
         m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w) /
        m_dv);
    tool.check_output(unorm, vnorm, wnorm);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.output_actions();
    // Only base level exists, which includes factor
    unorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u / m_dv);
    vnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v / m_dv);
    wnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);
}

TEST_F(FieldNormsTest, levelmask_on_int_refinement)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_intfieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_int_field(m_ifname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2 * m_ifref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5_rt;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.output_actions();
    // Only highest level will be counted
    // sqrt(number of cells * cell volume * cell value * cell value)
    amrex::Real unorm =
        std::sqrt((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_u * m_u / m_dv);
    amrex::Real vnorm =
        std::sqrt((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_v * m_v / m_dv);
    amrex::Real wnorm =
        std::sqrt((m_ncell0 * 8.0_rt) * (m_cv0 / 8.0_rt) * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.output_actions();
    // Only base level exists, which includes factor
    unorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u / m_dv);
    vnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v / m_dv);
    wnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);
}

TEST_F(FieldNormsTest, levelmask_not_cc)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& pressure = repo.declare_nd_field("pressure", 1, 1);
    auto& u_mac = repo.declare_xf_field("u_mac", 1, 1);
    auto& v_mac = repo.declare_yf_field("v_mac", 1, 1);
    auto& w_mac = repo.declare_zf_field("w_mac", 1, 1);
    auto& flag = repo.declare_field(m_fname, 1, 1);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2 * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize pressure distribution and access sim
    const amrex::Real fval = 10000.3_rt;
    pressure.setVal(fval);
    u_mac.setVal(fval);
    v_mac.setVal(fval);
    w_mac.setVal(fval);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up output (plot) variables
    io_mgr.register_output_var("pressure");
    io_mgr.register_output_var("u_mac");
    io_mgr.register_output_var("v_mac");
    io_mgr.register_output_var("w_mac");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.output_actions();
    // Only highest level will be counted
    tool.check_output4(fval);

    // Turn off refinements
    flag.setVal(0);
    rmesh.remesh();

    // Check again
    tool.output_actions();
    tool.check_output4(fval);
}

TEST_F(FieldNormsTest, norm_types)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_intfieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);

    // Initialize velocity distribution and access sim
    const amrex::Real factor = 0.1_rt;
    init_velocity(factor, velocity, m_u, m_v, m_w);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result
    FieldNormsImpl tool_l2(rsim, "fieldnorm");
    amrex::ParmParse pp("fieldnorm");
    pp.add("norm_type", (std::string) "2");
    tool_l2.initialize();
    tool_l2.output_actions();

    const amrex::Real l2_factor = std::sqrt(
        0.5_rt * ((1.0_rt - factor) * (1.0_rt - factor) +
                  (1.0_rt + factor) * (1.0_rt + factor)));
    amrex::Real unorm = m_u * l2_factor;
    amrex::Real vnorm = m_v * l2_factor;
    amrex::Real wnorm = m_w * l2_factor;
    tool_l2.check_output(unorm, vnorm, wnorm);

    FieldNormsImpl tool_l1(rsim, "fieldnorm");
    pp.add("norm_type", (std::string) "1");
    tool_l1.initialize();
    tool_l1.output_actions();

    const amrex::Real l1_factor =
        0.5_rt * ((1.0_rt - factor) + (1.0_rt + factor));
    unorm = m_u * l1_factor;
    vnorm = m_v * l1_factor;
    wnorm = m_w * l1_factor;
    tool_l1.check_output(unorm, vnorm, wnorm);

    FieldNormsImpl tool_linf(rsim, "fieldnorm");
    pp.add("norm_type", (std::string) "infinity");
    tool_linf.initialize();
    tool_linf.output_actions();

    const amrex::Real linf_factor = (1.0_rt + factor);
    unorm = m_u * linf_factor;
    vnorm = m_v * linf_factor;
    wnorm = m_w * linf_factor;
    tool_linf.check_output(unorm, vnorm, wnorm);
}

TEST_F(FieldNormsTest, norm_vector_magnitude)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_intfieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0_rt);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);

    // Initialize velocity distribution and access sim
    const amrex::Real factor = 0.1_rt;
    init_velocity(factor, velocity, m_u, m_v, m_w);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result
    FieldNormsImpl tool_l2(rsim, "fieldnorm");
    amrex::ParmParse pp("fieldnorm");
    pp.add("norm_type", (std::string) "2");
    pp.add("use_vector_magnitude", true);
    tool_l2.initialize();
    tool_l2.output_actions();

    const amrex::Real l2_factor = std::sqrt(
        0.5_rt * ((1.0_rt - factor) * (1.0_rt - factor) +
                  (1.0_rt + factor) * (1.0_rt + factor)));
    const amrex::Real vmag = std::sqrt((m_u * m_u) + (m_v * m_v) + (m_w * m_w));
    amrex::Real vmag_norm = vmag * l2_factor;
    tool_l2.check_output(vmag_norm);

    FieldNormsImpl tool_l1(rsim, "fieldnorm");
    pp.add("norm_type", (std::string) "1");
    tool_l1.initialize();
    tool_l1.output_actions();

    const amrex::Real l1_factor =
        0.5_rt * ((1.0_rt - factor) + (1.0_rt + factor));
    vmag_norm = vmag * l1_factor;
    tool_l1.check_output(vmag_norm);

    FieldNormsImpl tool_linf(rsim, "fieldnorm");
    pp.add("norm_type", (std::string) "infinity");
    tool_linf.initialize();
    tool_linf.output_actions();

    const amrex::Real linf_factor = (1.0_rt + factor);
    vmag_norm = vmag * linf_factor;
    tool_linf.check_output(vmag_norm);
}

} // namespace amr_wind_tests
