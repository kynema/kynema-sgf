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
void test_modes_distribution(bool ocean, bool cubic_time)
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
        data.cubic_time = cubic_time;
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
                        if (cubic_time) {
                            const double t = base + row * 108 + 3 * 18 + 2 * x;
                            const double fst = t + 2 * 18;
                            EXPECT_EQ(
                                data.c_mT[index],
                                std::complex<amrex::Real>(t, t + 1));
                            EXPECT_EQ(
                                data.c_mFST[index],
                                std::complex<amrex::Real>(fst, fst + 1));
                        }
                    }
                }
                if (cubic_time) {
                    EXPECT_EQ(data.c_mXt.size(), data.c_mT.size());
                    // Mean mode has no velocity or velocity rate
                    EXPECT_EQ(data.c_mXt[0], std::complex<amrex::Real>{});
                    EXPECT_EQ(data.c_mZt[0], std::complex<amrex::Real>{});
                } else {
                    EXPECT_TRUE(data.c_mT.empty());
                    EXPECT_TRUE(data.c_mXt.empty());
                }
            } else {
                const std::array arrays{
                    &data.r_mX, &data.r_mY, &data.r_mZ, &data.r_mFS};
                for (int row = 0; row < 2; ++row) {
                    for (int x = 0; x < 16; ++x) {
                        const int index = row * 16 + x;
                        for (int var = 0; var < 4; ++var) {
                            EXPECT_EQ(
                                (*arrays[var])[index],
                                base + row * 128 + (var == 3 ? 4 : var) * 16 +
                                    x);
                        }
                    }
                    for (int x = 0; x < 4; ++x) {
                        EXPECT_EQ(
                            data.r_mAdd[row * 4 + x],
                            base + row * 128 + 96 + x);
                        if (cubic_time) {
                            EXPECT_EQ(
                                data.r_mAddT[row * 4 + x],
                                base + row * 128 + 112 + x);
                        }
                    }
                    if (cubic_time) {
                        for (int x = 0; x < 16; ++x) {
                            const int index = row * 16 + x;
                            EXPECT_EQ(
                                data.r_mT[index], base + row * 128 + 48 + x);
                            EXPECT_EQ(
                                data.r_mFST[index], base + row * 128 + 80 + x);
                        }
                    }
                }
                if (cubic_time) {
                    EXPECT_EQ(data.r_mXt.size(), data.r_mT.size());
                    EXPECT_EQ(data.r_mXt[0], 0.0_rt);
                } else {
                    EXPECT_TRUE(data.r_mT.empty());
                    EXPECT_TRUE(data.r_mXt.empty());
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
    if (amrex::ParallelDescriptor::IOProcessor()) {
        EXPECT_EQ(std::remove(filename.c_str()), 0);
    }
}
} // namespace

TEST_F(OceanWavesW2ATest, modes_io_ocean)
{
    test_modes_distribution(true, false);
    test_modes_distribution(true, true);
}
TEST_F(OceanWavesW2ATest, modes_io_nwt)
{
    test_modes_distribution(false, false);
    test_modes_distribution(false, true);
}

class OceanWavesW2ACurrentTest : public OceanWavesW2ATest
{
protected:
    void check_current(
        bool init_waves,
        bool cubic_time,
        amrex::Real current,
        bool beach,
        bool waves = true)
    {
        populate_parameters();
        amrex::ParmParse geometry("geometry");
        geometry.addarr(
            "prob_lo", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, -4.0_rt});
        geometry.addarr(
            "prob_hi", amrex::Vector<amrex::Real>{8.0_rt, 8.0_rt, 4.0_rt});
        initialize_mesh();

        const auto filename =
            (std::filesystem::temp_directory_path() /
             ("w2a-current-" + amrex::UniqueString() + ".dat"))
                .string();
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream file(filename, std::ios::binary);
            const std::array<double, 10> header{16, 2, 1,    2, 8,
                                                8,  4, 9.81, 1, 1};
            for (int record = 0; record < 4; ++record) {
                for (int slot = 0; slot < 12 * 9 * 2; ++slot) {
                    double value = 0.0;
                    if (!record && slot < static_cast<int>(header.size())) {
                        value = header[slot];
                    } else if (record && waves && slot == 2) {
                        value = 0.2;
                    } else if (record && waves && slot == 4 * 18 + 2) {
                        value = 0.1;
                    }
                    file << std::uppercase << std::scientific
                         << std::setprecision(10) << std::setw(17) << value
                         << ' ';
                }
            }
            AMREX_ALWAYS_ASSERT(file.good());
        }

        amrex::ParmParse pp("OceanWaves.current_test");
        pp.add("HOS_modes_filename", filename);
        pp.add("current", current);
        pp.add("initialize_wave_field", init_waves);
        pp.add("number_interp_points_in_z", 10);
        pp.add("interp_spacing_at_surface", 0.5_rt);
        pp.add("number_interp_above_surface", 4);
        pp.add("relax_zone_gen_length", 1.0_rt);
        pp.add(
            beach ? "numerical_beach_length" : "relax_zone_out_length", 2.0_rt);
        pp.add("fft_reconstruction", std::string("full"));
        pp.add(
            "time_interpolation", std::string(cubic_time ? "cubic" : "linear"));
        auto& repo = sim().repo();
        auto& velocity = repo.declare_cc_field("velocity", 3, 3);
        auto& target = repo.declare_cc_field("ow_velocity", 3, 3);
        repo.declare_cc_field("levelset", 1, 3);
        repo.declare_cc_field("ow_levelset", 1, 3);
        velocity.setVal(amrex::Vector<amrex::Real>{8.0_rt, 9.0_rt, 10.0_rt}, 3);

        {
            using namespace kynema_sgf::ocean_waves;
            W2AWaves::DataType data(sim(), "current_test", 0);
            ::kynema_sgf::utils::MultiParser inputs(
                "OceanWaves.W2AWaves", "OceanWaves.current_test");
            ops::ReadInputsOp<W2AWaves>()(data, inputs);
            ops::InitDataOp<W2AWaves>()(data, 0, mesh().Geom(0), true);

            auto& cache = repo.get_field("w2a_velocity");
            auto& old_cache = cache.state(kynema_sgf::FieldState::Old);
            const amrex::Real tol = 1.0e-10_rt;
            const auto wave = utils::field_probe(old_cache, 0, 3, 3, 1);
            if (waves) {
                EXPECT_GT(std::abs(wave), 0.01_rt);
            } else {
                EXPECT_NEAR(wave, 0.0_rt, tol);
                ops::UpdateTargetFieldsOp<W2AWaves>()(data, 0.25_rt);
                EXPECT_NEAR(
                    utils::field_probe(target, 0, 3, 3, 1), current, tol);
            }
            EXPECT_NEAR(
                utils::field_probe(velocity, 0, 3, 3, 1),
                current + (init_waves ? wave : 0.0_rt), tol);
            EXPECT_NEAR(
                utils::field_probe(target, 0, 3, 3, 1), current + wave, tol);
            EXPECT_NEAR(utils::field_probe(velocity, 0, 3, 3, 7), 8.0_rt, tol);
            EXPECT_NEAR(
                utils::field_probe(velocity, 0, 3, 3, 7, 1), 9.0_rt, tol);
            EXPECT_NEAR(
                utils::field_probe(velocity, 0, 3, 3, 1, 1), 0.0_rt, tol);
            EXPECT_NEAR(utils::field_min(cache), 0.0_rt, tol);
            EXPECT_NEAR(utils::field_max(cache), 0.0_rt, tol);

            for (const amrex::Real time : {1.25_rt, 1.5_rt, 1.5_rt, 0.5_rt}) {
                if (time < data.meta().t_last) {
                    // Reconstruct both bracketing records, as after a restart.
                    data.meta().t_last = -1.0_rt;
                }
                ops::UpdateTargetFieldsOp<W2AWaves>()(data, time);
                EXPECT_NEAR(
                    utils::field_probe(target, 0, 3, 3, 1), current + wave,
                    tol);
                EXPECT_NEAR(utils::field_probe(cache, 0, 3, 3, 1), wave, tol);
                EXPECT_NEAR(
                    utils::field_probe(old_cache, 0, 3, 3, 1), wave, tol);
                EXPECT_NEAR(utils::field_min(target, 1), 0.0_rt, tol);
                EXPECT_NEAR(utils::field_max(target, 2), 0.0_rt, tol);
                if (beach) {
                    const auto before = utils::field_probe(target, 0, 7, 3, 1);
                    relaxation_zones::modify_target_fields_for_beach(
                        sim(), data.meta());
                    const auto after = utils::field_probe(target, 0, 7, 3, 1);
                    if (waves) {
                        EXPECT_LT(
                            std::abs(after - current),
                            std::abs(before - current));
                    } else {
                        EXPECT_NEAR(after, current, tol);
                    }
                }
            }
        }
        if (amrex::ParallelDescriptor::IOProcessor()) {
            EXPECT_EQ(std::remove(filename.c_str()), 0);
        }
    }
};

TEST_F(OceanWavesW2ACurrentTest, following_current_with_waves)
{
    check_current(true, false, 1.25_rt, false);
}

TEST_F(OceanWavesW2ACurrentTest, opposing_current_without_waves)
{
    check_current(false, true, -0.75_rt, true);
}

TEST_F(OceanWavesW2ACurrentTest, zero_current_without_waves)
{
    check_current(false, false, 0.0_rt, true);
}

TEST_F(OceanWavesW2ACurrentTest, current_only)
{
    check_current(false, true, 1.0_rt, true, false);
}
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
    auto& w2a_lvs = sim().repo().declare_cc_field("w2a_levelset", 1, 3, 2);
    auto& w2a_vel = sim().repo().declare_cc_field("w2a_velocity", 3, 3, 2);
    auto& w2a_lvs_old = w2a_lvs.state(kynema_sgf::FieldState::Old);
    auto& w2a_vel_old = w2a_vel.state(kynema_sgf::FieldState::Old);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;

    // Linear: old data at t0 and new data at t1
    w2a_lvs_old.setVal(1.0_rt);
    w2a_lvs.setVal(0.0_rt);
    w2a_vel_old.setVal(-2.0_rt);
    w2a_vel.setVal(2.0_rt);
    time_interpolate_wave_fields(
        lvs, vel, w2a_lvs_old, w2a_vel_old, w2a_lvs, w2a_vel, nullptr, nullptr,
        nullptr, nullptr, 1.0_rt, 2.0_rt, 1.5_rt);
    EXPECT_NEAR(utils::field_max(lvs), 0.5_rt, tol);
    EXPECT_NEAR(utils::field_min(lvs), 0.5_rt, tol);
    EXPECT_NEAR(utils::field_max(vel), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_min(vel), 0.0_rt, tol);
    // Outside of the interval, data is held rather than extrapolated
    time_interpolate_wave_fields(
        lvs, vel, w2a_lvs_old, w2a_vel_old, w2a_lvs, w2a_vel, nullptr, nullptr,
        nullptr, nullptr, 1.0_rt, 2.0_rt, 2.5_rt);
    EXPECT_NEAR(utils::field_max(lvs), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_max(vel), 2.0_rt, tol);

    // Cubic Hermite is exact for data that is cubic in time
    auto& dlvs = sim().repo().declare_cc_field("w2a_dlevelset_dt", 1, 3, 2);
    auto& dvel = sim().repo().declare_cc_field("w2a_dvelocity_dt", 3, 3, 2);
    auto& dlvs_old = dlvs.state(kynema_sgf::FieldState::Old);
    auto& dvel_old = dvel.state(kynema_sgf::FieldState::Old);
    const auto f = [](amrex::Real t) {
        return 0.7_rt * t * t * t - 2.0_rt * t * t + 0.5_rt * t + 1.0_rt;
    };
    const auto df = [](amrex::Real t) {
        return 2.1_rt * t * t - 4.0_rt * t + 0.5_rt;
    };
    const amrex::Real t0 = 0.25_rt;
    const amrex::Real t1 = 1.75_rt;
    w2a_lvs_old.setVal(f(t0));
    w2a_lvs.setVal(f(t1));
    dlvs_old.setVal(df(t0));
    dlvs.setVal(df(t1));
    w2a_vel_old.setVal(-f(t0));
    w2a_vel.setVal(-f(t1));
    dvel_old.setVal(-df(t0));
    dvel.setVal(-df(t1));
    for (const amrex::Real t : {t0, 0.6_rt, 1.0_rt, 1.3_rt, t1}) {
        time_interpolate_wave_fields(
            lvs, vel, w2a_lvs_old, w2a_vel_old, w2a_lvs, w2a_vel, &dlvs_old,
            &dvel_old, &dlvs, &dvel, t0, t1, t);
        EXPECT_NEAR(utils::field_max(lvs), f(t), tol);
        EXPECT_NEAR(utils::field_min(lvs), f(t), tol);
        EXPECT_NEAR(utils::field_max(vel), -f(t), tol);
        EXPECT_NEAR(utils::field_min(vel), -f(t), tol);
    }

    // Hermite weights reduce to the endpoint values and partition unity
    for (const amrex::Real t : {t0, 0.9_rt, t1}) {
        const auto w = time_interp_weights(t0, t1, t, true);
        EXPECT_NEAR(w[0] + w[1], 1.0_rt, tol);
        const auto l = time_interp_weights(t0, t1, t, false);
        EXPECT_NEAR(l[0] + l[1], 1.0_rt, tol);
        EXPECT_EQ(l[2], 0.0_rt);
        EXPECT_EQ(l[3], 0.0_rt);
    }
    const auto w0 = time_interp_weights(t0, t1, t0, true);
    EXPECT_EQ(w0[0], 1.0_rt);
    EXPECT_EQ(w0[2], 0.0_rt);
    EXPECT_EQ(w0[3], 0.0_rt);
}

} // namespace kynema_sgf_tests
