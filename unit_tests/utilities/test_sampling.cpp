#include <numbers>
#include "ks_test_utils/MeshTest.H"
#include "src/utilities/sampling/Sampling.H"
#include "src/utilities/sampling/SamplingContainer.H"
#include "src/utilities/sampling/ProbeSampler.H"
#include "src/utilities/sampling/PlaneSampler.H"
#include "src/utilities/sampling/VolumeSampler.H"
#include "src/utilities/sampling/DTUSpinnerSampler.H"
#include "src/utilities/sampling/RadarSampler.H"
#include "src/utilities/sampling/SamplingUtils.H"
#include "AMReX_Vector.H"
#include "src/core/vs/vector_space.H"
#include "src/utilities/tensor_ops.H"
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

void init_int_field(kynema_sgf::IntField& fld)
{
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(), ncomp,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                farrs[nbx](i, j, k, n) = i + j + k;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void write_probe_sampler_file(const std::string& fname)
{
    std::ofstream os(fname);
    // Total number of points
    os << "3\n";
    // Coordinates
    os << "0.0\t0.0\t0.0\n";
    os << "60.0\t2.0\t3.0\n";
    os << "100.0\t8.0\t5.0\n";
}

class SamplingImpl : public kynema_sgf::sampling::Sampling
{
public:
    SamplingImpl(kynema_sgf::CFDSim& sim, const std::string& label)
        : kynema_sgf::sampling::Sampling(sim, label)
    {}

    bool write_flag{false};

protected:
    void prepare_netcdf_file() override {}
    void process_output() override
    {
        // Test buffer populate for GPU runs
        std::vector<amrex::Real> buf(
            num_total_particles() * var_names().size(), 0.0_rt);
        sampling_container().populate_buffer(buf);

        write_flag = true;
    }
};

} // namespace

class SamplingTest : public MeshTest
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

namespace {

void test_scontainer_impl(
    kynema_sgf::sampling::SamplingContainer::ParticleVector& pvec,
    const int npts)
{
    using IIx = kynema_sgf::sampling::IIx;
    using PType = kynema_sgf::sampling::SamplingContainer::ParticleType;
    // Create a line probe
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> begin{
        {66.0_rt, 66.0_rt, 1.0_rt}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> end{
        {66.0_rt, 66.0_rt, 127.0_rt}};

    auto* pstruct = pvec.data();
    const int id_start = static_cast<int>(PType::NextID());
    const int proc_id = amrex::ParallelDescriptor::MyProc();

    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE(const int ip) {
        auto& pp = pstruct[ip];
        pp.id() = id_start + ip;
        pp.cpu() = proc_id;
        pp.pos(0) = begin[0];
        pp.pos(1) = begin[1];

        const int dl = static_cast<int>(
            (end[2] - begin[2]) / static_cast<amrex::Real>(npts - 1));
        const amrex::Real z = (ip + 0.5_rt) * dl;
        pp.pos(2) = z;

        pp.idata(IIx::uid) = static_cast<int>(pp.id());
        pp.idata(IIx::sid) = 0;
        pp.idata(IIx::nid) = ip;
    });
    amrex::Gpu::streamSynchronize();
}

} // namespace

TEST_F(SamplingTest, scontainer)
{
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        GTEST_SKIP();
    }

    using ParIter = kynema_sgf::sampling::SamplingContainer::ParIterType;
    const int ncomp = 5;
    initialize_mesh();

    kynema_sgf::sampling::SamplingContainer sc(mesh());
    // Set up runtime components and prime the particle container
    sc.setup_container(ncomp);

    const int npts = 64;
    const int lev = 0;
    const int gid = 0;
    const int tid = 0;
    auto& ptile = sc.GetParticles(lev)[std::make_pair(gid, tid)];
    ptile.resize(npts);
    test_scontainer_impl(ptile.GetArrayOfStructs()(), npts);
    sc.Redistribute();

#ifndef AMREX_USE_GPU
    // TODO: Check why this causes errors in non-managed CUDA run
    ASSERT_EQ(npts, sc.NumberOfParticlesAtLevel(lev));
#endif
    int counter = 0;
    int total_particles = 0;
    for (ParIter pti(sc, lev); pti.isValid(); ++pti) {
        ++counter;
        total_particles += pti.numParticles();
    }
    ASSERT_EQ(total_particles, npts);
    ASSERT_EQ(counter, 4);
}

TEST_F(SamplingTest, sampling)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vel = repo.declare_field("velocity", 3, 2);
    auto& pres = repo.declare_nd_field("pressure", 1, 2);
    auto& rho = repo.declare_field("density", 1, 2);
    auto& idxsum = repo.declare_int_field("idxsum", 1, 2);
    init_field(vel);
    init_field(pres);
    init_field(rho);
    init_int_field(idxsum);

    {
        amrex::ParmParse pp("sampling");
        pp.add("output_interval", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"line1"});
        pp.addarr(
            "fields",
            amrex::Vector<std::string>{"density", "pressure", "velocity"});
        pp.addarr(
            "derived_fields", amrex::Vector<std::string>{"mag_vorticity"});
        pp.addarr("int_fields", amrex::Vector<std::string>{"idxsum"});
    }
    {
        amrex::ParmParse pp("sampling.line1");
        pp.add("type", std::string("LineSampler"));
        pp.add("num_points", 16);
        pp.addarr(
            "start", amrex::Vector<amrex::Real>{66.0_rt, 66.0_rt, 1.0_rt});
        pp.addarr(
            "end", amrex::Vector<amrex::Real>{66.0_rt, 66.0_rt, 127.0_rt});
    }

    SamplingImpl probes(sim(), "sampling");
    probes.initialize();
    probes.output_actions();

    EXPECT_TRUE(probes.write_flag);
}

TEST_F(SamplingTest, sampling_timing)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vel = repo.declare_field("velocity", 3, 2);
    auto& pres = repo.declare_nd_field("pressure", 1, 2);
    auto& rho = repo.declare_field("density", 1, 2);
    auto& idxsum = repo.declare_int_field("idxsum", 1, 2);
    init_field(vel);
    init_field(pres);
    init_field(rho);
    init_int_field(idxsum);

    {
        amrex::ParmParse pp("sampling");
        pp.add("output_interval", 1);
        pp.add("output_delay", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"line1"});
        pp.addarr(
            "fields",
            amrex::Vector<std::string>{"density", "pressure", "velocity"});
        pp.addarr(
            "derived_fields", amrex::Vector<std::string>{"mag_vorticity"});
        pp.addarr("int_fields", amrex::Vector<std::string>{"idxsum"});
    }
    {
        amrex::ParmParse pp("sampling.line1");
        pp.add("type", std::string("LineSampler"));
        pp.add("num_points", 16);
        pp.addarr(
            "start", amrex::Vector<amrex::Real>{66.0_rt, 66.0_rt, 1.0_rt});
        pp.addarr(
            "end", amrex::Vector<amrex::Real>{66.0_rt, 66.0_rt, 127.0_rt});
    }

    SamplingImpl probes(sim(), "sampling");
    probes.initialize();
    if (probes.do_output_now(
            sim().time().time_index(), sim().time().start_time_index(),
            sim().time().new_time(), sim().time().start_time(),
            sim().time().delta_t(), 1.0_rt)) {
        probes.output_actions();
    }
    EXPECT_FALSE(probes.write_flag);
    sim().time().new_timestep();
    if (probes.do_output_now(
            sim().time().time_index(), sim().time().start_time_index(),
            sim().time().new_time(), sim().time().start_time(),
            sim().time().delta_t(), 1.0_rt)) {
        probes.output_actions();
    }
    EXPECT_TRUE(probes.write_flag);
}

TEST_F(SamplingTest, probe_sampler)
{
    initialize_mesh();

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    std::string fname = "probes.txt";
    // Write file
    write_probe_sampler_file(fname);

    {
        amrex::ParmParse pp("cloud");
        pp.add("probe_location_file", fname);
        pp.addarr("offsets", amrex::Vector<amrex::Real>{1.0_rt, 2.5_rt});
        pp.addarr(
            "offset_vector",
            amrex::Vector<amrex::Real>{0.2_rt, 0.5_rt, 1.0_rt});
    }

    kynema_sgf::sampling::ProbeSampler cloud(sim());
    cloud.initialize("cloud");
    kynema_sgf::sampling::SampleLocType sample_locs;
    cloud.sampling_locations(sample_locs);

    ASSERT_EQ(sample_locs.locations().size(), 3 * 2);
    const amrex::Vector<amrex::Real> xprobe_golds{0.2_rt, 60.2_rt, 100.2_rt,
                                                  0.5_rt, 60.5_rt, 100.5_rt};
    const amrex::Vector<amrex::Real> yprobe_golds{0.5_rt,  2.5_rt,  8.5_rt,
                                                  1.25_rt, 3.25_rt, 9.25_rt};
    const amrex::Vector<amrex::Real> zprobe_golds{1.0_rt, 4.0_rt, 6.0_rt,
                                                  2.5_rt, 5.5_rt, 7.5_rt};
    const auto& locs = sample_locs.locations();
    for (int n = 0; n < locs.size(); ++n) {
        EXPECT_NEAR(locs[n][0], xprobe_golds[n], tol);
        EXPECT_NEAR(locs[n][1], yprobe_golds[n], tol);
        EXPECT_NEAR(locs[n][2], zprobe_golds[n], tol);
    }

    // Remove file
    const char* fname_char = fname.c_str();
    {
        std::ifstream f(fname_char);
        if (f.good()) {
            remove(fname_char);
        }
        // Check that file is removed
        std::ifstream ff(fname_char);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(SamplingTest, plane_sampler)
{
    initialize_mesh();

    {
        amrex::ParmParse pp("plane");
        pp.addarr("axis1", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt, 0.0_rt});
        pp.addarr("axis2", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 1.0_rt});
        pp.addarr("origin", amrex::Vector<amrex::Real>{1.0_rt, 1.0_rt, 1.0_rt});
        pp.addarr("num_points", amrex::Vector<int>{3, 3});
        pp.addarr("offsets", amrex::Vector<amrex::Real>{2.0_rt, 10.0_rt});
        pp.addarr(
            "offset_vector",
            amrex::Vector<amrex::Real>{1.0_rt, 0.0_rt, 0.0_rt});
    }

    kynema_sgf::sampling::PlaneSampler plane(sim());
    plane.initialize("plane");
    kynema_sgf::sampling::SampleLocType sample_locs;
    plane.sampling_locations(sample_locs);

    ASSERT_EQ(sample_locs.locations().size(), 3 * 3 * 2);
}

TEST_F(SamplingTest, volume_sampler)
{
    initialize_mesh();
    amrex::ParmParse pp("volume");
    pp.addarr("hi", amrex::Vector<amrex::Real>{1.0_rt, 1.0_rt, 1.0_rt});
    pp.addarr("lo", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
    pp.addarr("num_points", amrex::Vector<int>{3, 5, 5});

    kynema_sgf::sampling::VolumeSampler volume(sim());
    volume.initialize("volume");
    kynema_sgf::sampling::SampleLocType sample_locs;
    volume.sampling_locations(sample_locs);

    ASSERT_EQ(sample_locs.locations().size(), 3 * 5 * 5);
}

TEST_F(SamplingTest, spinner_sampler)
{
    initialize_mesh();
    amrex::ParmParse pp("spinner");

    pp.add("mode", std::string("fixed"));
    pp.add("turbine", std::string("WTG01"));
    pp.add("hub_debug", false);
    pp.add("inner_prism_theta0", 90.0_rt);
    pp.add("inner_prism_rotrate", 3.5_rt);
    pp.add("inner_prism_azimuth", 15.2_rt);
    pp.add("outer_prism_theta0", 90.0_rt);
    pp.add("outer_prism_rotrate", 6.5_rt);
    pp.add("outer_prism_azimuth", 15.2_rt);
    pp.addarr(
        "lidar_center",
        amrex::Vector<amrex::Real>{630.0_rt, 192.0_rt, 120.0_rt});
    pp.add("scan_time", 2.0_rt);
    pp.add("num_samples", 984);
    pp.add("beam_length", 270.0_rt);
    pp.add("beam_points", 432);
    pp.add("fixed_yaw", 0);
    pp.add("fixed_roll", 0);
    pp.add("fixed_tilt", 0);

    kynema_sgf::sampling::DTUSpinnerSampler spinner(sim());
    spinner.initialize("spinner");
    kynema_sgf::sampling::SampleLocType sample_locs;
    spinner.sampling_locations(sample_locs);

    ASSERT_EQ(sample_locs.locations().size(), 21600);
}

TEST_F(SamplingTest, radar_sampler)
{
    initialize_mesh();
    amrex::ParmParse pp("radar");

    pp.add("num_points", 512);
    pp.addarr("origin", amrex::Vector<amrex::Real>{1.0_rt, 1.0_rt, 1.0_rt});
    pp.add("sampling_frequency", 85.0_rt);
    pp.add("device_sampling_frequency", 30.0_rt);
    pp.add("radar_cone_angle", 0.25_rt);
    pp.add("radar_quadrature_type", std::string("truncated_normal_halfpower"));
    pp.add("radar_npts_azimuth", 5);
    pp.add("radar_beam_length", 100.0_rt);
    pp.add("angular_speed", 30.0_rt);
    pp.add("sweep_angle", 145.0_rt);
    pp.add("reset_time", 0.0_rt);
    pp.addarr(
        "elevation_angles",
        amrex::Vector<amrex::Real>{0.0_rt, 0.1_rt, 0.2_rt, 0.3_rt, 0.4_rt});
    pp.addarr(
        "axis",
        amrex::Vector<amrex::Real>{0.707106781_rt, 0.707106781_rt, 0.0_rt});
    pp.addarr(
        "vertical_unit_dir",
        amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 1.0_rt});
    pp.add("debug_print", false);

    kynema_sgf::sampling::RadarSampler radar(sim());
    radar.initialize("radar");
    kynema_sgf::sampling::SampleLocType sample_locs;
    radar.sampling_locations(sample_locs);

    ASSERT_EQ(sample_locs.locations().size(), 193536);
}

TEST_F(SamplingTest, sampling_utils)
{
    namespace vs = kynema_sgf::vs;
    amrex::Real toler = std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;
    vs::Vector unitx{1.0_rt, 0.0_rt, 0.0_rt};
    vs::Vector unity{0.0_rt, 1.0_rt, 0.0_rt};
    vs::Vector unitz{0.0_rt, 0.0_rt, 1.0_rt};
    vs::Vector nunitx{-1.0_rt, 0.0_rt, 0.0_rt};
    vs::Vector ffn{-0.70710678_rt, -0.70710678_rt, 0.0_rt};
    vs::Vector ffnr{-0.70710678_rt, 0.70710678_rt, 0.0_rt};
    vs::Vector result;
    vs::Vector angles{0.0_rt, 180.0_rt, 0.0_rt};

    result = kynema_sgf::sampling::sampling_utils::reflect(unity, ffn);
    EXPECT_NEAR(result[0], ffnr[0], toler);
    EXPECT_NEAR(result[1], ffnr[1], toler);
    EXPECT_NEAR(result[2], ffnr[2], toler);

    result = kynema_sgf::sampling::sampling_utils::rotation(angles, unitx);
    EXPECT_NEAR(result[0], nunitx[0], toler);
    EXPECT_NEAR(result[1], nunitx[1], toler);
    EXPECT_NEAR(result[2], nunitx[2], toler);

    result = kynema_sgf::sampling::sampling_utils::rotate_euler_vec(
        unity, -90.0_rt, unitx);
    EXPECT_NEAR(result[0], unitz[0], toler);
    EXPECT_NEAR(result[1], unitz[1], toler);
    EXPECT_NEAR(result[2], unitz[2], toler);
}

TEST_F(SamplingTest, quadrature)
{
    amrex::Real toler = std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    namespace vs = kynema_sgf::vs;
    namespace su = kynema_sgf::sampling::sampling_utils;
    int ntheta = 5;
    amrex::Real gammav = 0.25_rt * std::numbers::pi_v<amrex::Real> / 180.0_rt;
    std::vector<amrex::Real> weights = {
        0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt,
        0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt,
        0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt};
    vs::Vector sr{0.0_rt, 0.0_rt, 1.0_rt};
    std::vector<vs::Vector> rays = {sr, sr, sr, sr, sr, sr, sr, sr, sr, sr, sr,
                                    sr, sr, sr, sr, sr, sr, sr, sr, sr, sr};

    su::spherical_cap_truncated_normal(
        gammav, ntheta, su::NormalRule::HALFPOWER, rays, weights);

    EXPECT_NEAR(weights[0], 1.1826123083219489e-05_rt, toler);
    EXPECT_NEAR(weights[10], 3.004263660003298e-06_rt, toler);
    EXPECT_NEAR(weights[20], 6.6402168628164281e-07_rt, toler);
}

} // namespace kynema_sgf_tests
