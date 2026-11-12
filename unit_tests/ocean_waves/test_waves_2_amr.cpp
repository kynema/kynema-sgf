#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/ocean_waves/utils/wave_utils_K.H"
#include "src/ocean_waves/OceanWaves.H"
#include "src/ocean_waves/relaxation_zones/waves2amr_ops.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"
#include <filesystem>
#include <fstream>
#include <iomanip>

using namespace amrex::literals;

namespace kynema_sgf_tests {

class OceanWavesW2ATest : public MeshTest
{};

#ifdef KYNEMA_SGF_USE_W2A
namespace {
void test_modes_distribution(bool ocean)
{
    const auto filename = (std::filesystem::temp_directory_path() /
                           ("w2a-mpi-" + amrex::UniqueString() + ".dat"))
                              .string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream file(filename, std::ios::binary);
        const int width = ocean ? 9 : 16;
        const int slots = (ocean ? 12 : 8) * width * 2;
        const std::vector<double> header =
            ocean ? std::vector<double>{16, 2, .1, 1, 20, 10, 5, 9.81, 1, 1}
                  : std::vector<double>{16, 4, 2, .1, 1, 20, 10, 5, 0};
        for (int record = 0; record < 4; ++record) {
            for (int slot = 0; slot < slots; ++slot) {
                double value = record * 1000 + slot;
                if (!record && slot < static_cast<int>(header.size()))
                    value = header[slot];
                if (!ocean && slot % (8 * width) >= 6 * width &&
                    slot % width >= 4) {
                    file << std::string(18, ' ');
                } else {
                    file << std::uppercase << std::scientific
                         << std::setprecision(10) << std::setw(17) << value
                         << ' ';
                }
            }
        }
        AMREX_ALWAYS_ASSERT(file.good());
    }
    {
        kynema_sgf::ocean_waves::W2AWavesData data;
        data.is_ocean = ocean;
        // The I/O rank is deliberately inactive; only it may open the file.
        data.do_interp = !amrex::ParallelDescriptor::IOProcessor();
        if (ocean) {
            initialize_modes_reader(data.c_rmodes, filename, true);
            for (auto* modes :
                 {&data.c_mX, &data.c_mY, &data.c_mZ, &data.c_mFS})
                modes->resize(data.c_rmodes.get_vector_size());
        } else {
            initialize_modes_reader(data.r_rmodes, filename, false);
            for (auto* modes :
                 {&data.r_mX, &data.r_mY, &data.r_mZ, &data.r_mFS})
                modes->resize(data.r_rmodes.get_vector_size());
            data.r_mAdd.resize(data.r_rmodes.get_addl_vector_size());
        }
        const auto parses = [&] {
            return ocean ? data.c_rmodes.records_read()
                         : data.r_rmodes.records_read();
        };
        const auto expected_parses = [](int count) {
            return amrex::ParallelDescriptor::IOProcessor() ? count : 0;
        };
        const auto check = [&](int step) {
            EXPECT_EQ(data.loaded_modes_step, step);
            const double base = (step + 1) * 1000;
            if (ocean) {
                const std::array arrays{
                    &data.c_mX, &data.c_mY, &data.c_mZ, &data.c_mFS};
                for (int row = 0; row < 2; ++row) {
                    for (int x = 0; x < 9; ++x) {
                        const int index = row * 9 + x;
                        for (int var = 0; var < 4; ++var) {
                            const double real = base + row * 108 +
                                                (var == 3 ? 4 : var) * 18 +
                                                2 * x;
                            EXPECT_EQ(
                                (*arrays[var])[index],
                                std::complex<amrex::Real>(real, real + 1));
                        }
                    }
                }
            } else {
                const std::array arrays{
                    &data.r_mX, &data.r_mY, &data.r_mZ, &data.r_mFS};
                for (int row = 0; row < 2; ++row) {
                    for (int x = 0; x < 16; ++x) {
                        const int index = row * 16 + x;
                        for (int var = 0; var < 4; ++var)
                            EXPECT_EQ(
                                (*arrays[var])[index],
                                base + row * 128 + (var == 3 ? 4 : var) * 16 +
                                    x);
                    }
                    for (int x = 0; x < 4; ++x)
                        EXPECT_EQ(
                            data.r_mAdd[row * 4 + x],
                            base + row * 128 + 96 + x);
                }
            }
        };
        data.ntime = 2;
        read_wave_modes(data);
        check(2);
        EXPECT_EQ(parses(), expected_parses(1));
        read_wave_modes(data);
        check(2);
        EXPECT_EQ(parses(), expected_parses(1));
        // Backward/restart-style request, then EOF looping back to cached data.
        data.ntime = 0;
        read_wave_modes(data);
        check(0);
        data.ntime = 8;
        read_wave_modes(data);
        check(0);
        EXPECT_EQ(data.n_offset, -8);
        EXPECT_EQ(parses(), expected_parses(2));
        data.ntime = 9;
        read_wave_modes(data);
        check(1);
        EXPECT_EQ(parses(), expected_parses(3));
    }
    if (amrex::ParallelDescriptor::IOProcessor())
        EXPECT_EQ(std::remove(filename.c_str()), 0);
}
} // namespace

TEST_F(OceanWavesW2ATest, modes_io_ocean) { test_modes_distribution(true); }
TEST_F(OceanWavesW2ATest, modes_io_nwt) { test_modes_distribution(false); }
#endif

// This is to evaluate aspects of W2A methods that do not require wave data
TEST_F(OceanWavesW2ATest, time_advance)
{
    // Variables modified by routine
    int ntime = 0;
    bool read_flag = false;
    bool resize_flag = false;
    amrex::Real t = 0.0_rt;
    amrex::Real t_last = -1.0_rt;
    // Offset from wave data
    const amrex::Real t_winit = 0.0_rt;

    // Upon initialization with waves
    int new_ntime = 1;
    amrex::Real newtime = 0.001_rt;
    amrex::Real dt_modes = 0.1_rt;
    t_last = 0.0_rt; // Set within InitDataOp
    amrex::Real t_last_prior = t_last;
    int double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    // Time interpolation factor
    amrex::Real f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    // Should regrid and read, time interpolation should be ordinary. Read once
    EXPECT_EQ(double_data, 0);
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 0.1_rt);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, newtime / t,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // Upon initialization without waves, when t_last has not been modified
    read_flag = false;
    resize_flag = false;
    ntime = 0;
    newtime = 0.001_rt;
    t_last = -1.0_rt;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    EXPECT_EQ(double_data, 1);
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 0.1_rt);
    EXPECT_EQ(t_last, newtime);
    // Double data routine changes t_last_prior
    t_last_prior = (ntime - 1) * dt_modes;
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_NEAR(
        f_interp, newtime / t,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // Upon initialization as restart, when t_last has not been modified
    read_flag = false;
    resize_flag = false;
    t_last = -1.0_rt;
    ntime = 0;
    new_ntime = 40;
    newtime = 3.99_rt;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    EXPECT_EQ(double_data, 1);
    EXPECT_EQ(ntime, 40);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 40 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    // Double data routine changes t_last_prior
    t_last_prior = (ntime - 1) * dt_modes;
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_NEAR(
        f_interp, (3.99_rt - 3.9_rt) / dt_modes,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // Upon initialization as restart, when newtime happens to be at ntime
    read_flag = false;
    resize_flag = false;
    t_last = -1.0_rt;
    t_last_prior = t_last;
    ntime = 0;
    new_ntime = 40;
    newtime = 4.0_rt;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(double_data, 2);
    EXPECT_EQ(ntime, 40);
    EXPECT_TRUE(read_flag);
    EXPECT_TRUE(resize_flag);
    EXPECT_EQ(t, 40 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, 1.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // The following 3 tests are for during initial iterations (after the first)
    // -- Beginning from time = 0. --
    read_flag = false;
    resize_flag = false;
    ntime = 1;
    t = ntime * dt_modes;
    new_ntime = 1;
    newtime = 0.001_rt;
    t_last = 0.001_rt;
    t_last_prior = t_last;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    // ---- Should not regrid nor read, time interpolation should do nothing
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    // -- Beginning from time = 3.98_rt --
    read_flag = false;
    resize_flag = false;
    t_last = 3.99_rt;
    t_last_prior = t_last;
    ntime = 40;
    t = ntime * dt_modes;
    new_ntime = 40;
    newtime = 3.99_rt;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    // ---- Should not regrid nor read, time interpolation should do nothing
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    // -- Beginning from time = 3.99_rt --
    read_flag = false;
    resize_flag = false;
    t_last = 4.00_rt;
    t_last_prior = t_last;
    ntime = 40;
    t = ntime * dt_modes;
    new_ntime = 40;
    newtime = 4.00_rt;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) /
               (t - t_last_prior + std::numeric_limits<amrex::Real>::epsilon());
    // Should not regrid nor read, time interpolation should do nothing
    EXPECT_EQ(ntime, new_ntime);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // The following tests are for during ordinary timesteps
    // -- No resize, no reading, just interp
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    newtime = 3.99_rt;
    t = 40 * dt_modes;
    t_last = 3.98_rt;
    new_ntime = 40;
    t_last_prior = t_last;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(ntime, 40);
    EXPECT_FALSE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, 0.01_rt / 0.02_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    // -- No resize, yes reading
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    new_ntime = 41;
    newtime = 4.01_rt;
    t_last = 4.0_rt;
    t_last_prior = t_last;
    evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(ntime, 41);
    EXPECT_TRUE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t, 41 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    EXPECT_NEAR(
        f_interp, 0.01_rt / (4.1_rt - 4.0_rt),
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    // -- No resize, yes reading, big timestep means amrex::Real read
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    new_ntime = 42;
    newtime = 4.15_rt;
    t_last = 4.0_rt;
    // built-in change to t_last_prior
    t_last_prior = (new_ntime - 1) * dt_modes;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(double_data, 1);
    EXPECT_EQ(ntime, 42);
    EXPECT_TRUE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t, 42 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    // interp is from 4.1_rt (replaced) to 4.15_rt
    EXPECT_NEAR(
        f_interp, (4.15_rt - 4.1_rt) / (4.2_rt - 4.1_rt),
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    // -- No resize, yes reading, big timestep leads to convenient single read
    resize_flag = false;
    read_flag = false;
    ntime = 40;
    new_ntime = 42;
    newtime = 4.2_rt;
    t_last = 4.0_rt;
    t_last_prior = t_last;
    double_data = evaluate_read_resize(
        ntime, read_flag, resize_flag, t, t_last, new_ntime, t_winit, dt_modes,
        newtime);
    f_interp = (newtime - t_last_prior) / (t - t_last_prior);
    EXPECT_EQ(double_data, 2);
    EXPECT_EQ(ntime, 42);
    EXPECT_TRUE(read_flag);
    EXPECT_FALSE(resize_flag);
    EXPECT_EQ(t, 42 * dt_modes);
    EXPECT_EQ(t_last, newtime);
    // interp is from 4.0_rt to 4.2_rt, but it leads to 1
    EXPECT_NEAR(
        f_interp, (4.2_rt - 4.0_rt) / (4.2_rt - 4.0_rt),
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
}

TEST_F(OceanWavesW2ATest, time_reset)
{
    // Testing the arithmetic for resetting data back to the beginning
    int nstop = 10;
    int n0 = 1;

    int ntime = 11;
    int noff = update_offset_timestep(ntime, n0);
    EXPECT_EQ(ntime + noff, 1);

    // As if continued through another set of data
    ntime = nstop + ntime;
    noff = update_offset_timestep(ntime, n0);
    EXPECT_EQ(ntime + noff, 1);
}

TEST_F(OceanWavesW2ATest, time_interp)
{
    populate_parameters();
    initialize_mesh();
    // Set up fields and initial values
    auto& lvs = sim().repo().declare_cc_field("ow_levelset", 1, 3, 1);
    auto& vel = sim().repo().declare_cc_field("ow_velocity", 3, 3, 1);
    auto& w2a_lvs = sim().repo().declare_cc_field("w2a_levelset", 1, 3, 1);
    auto& w2a_vel = sim().repo().declare_cc_field("w2a_velocity", 3, 3, 1);
    lvs.setVal(0.0_rt);
    amrex::Real t_last = 0.0_rt;

    // First time from "file"
    w2a_lvs.setVal(1.0_rt);
    amrex::Real t_w2a = 1.0_rt;

    // Interpolate once in time
    amrex::Real t_sim = 0.5_rt;
    time_interpolate_wave_fields(
        lvs, vel, w2a_lvs, w2a_vel, t_last, t_sim, t_w2a);

    // Record time that data has been interpolated to
    t_last = t_sim;

    // Check value of levelset
    const amrex::Real fmax1 = utils::field_max(lvs);
    const amrex::Real fmin1 = utils::field_min(lvs);
    EXPECT_NEAR(
        fmax1, 0.5_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_NEAR(
        fmin1, 0.5_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);

    // Step forward in time, after current w2a data
    t_sim += 1.;

    if (t_sim > t_w2a) {
        // Copy current w2a data and update old time
        kynema_sgf::field_ops::copy(lvs, w2a_lvs, 0, 0, 1, lvs.num_grow());
        t_last = t_w2a;
        // Update w2a values (as in a file read)
        w2a_lvs.setVal(0.0_rt);
        t_w2a += 1.0_rt;
    }

    // Interpolate again
    time_interpolate_wave_fields(
        lvs, vel, w2a_lvs, w2a_vel, t_last, t_sim, t_w2a);

    // Check value of levelset
    const amrex::Real fmax2 = utils::field_max(lvs);
    const amrex::Real fmin2 = utils::field_min(lvs);
    EXPECT_NEAR(
        fmax2, 0.5_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_NEAR(
        fmin2, 0.5_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
}

} // namespace kynema_sgf_tests
