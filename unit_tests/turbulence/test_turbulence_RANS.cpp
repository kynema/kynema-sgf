#include <cmath>
#include <limits>
#include <string>

#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/turbulence/TurbulenceModel.H"
#include "ks_test_utils/test_utils.H"
#include "src/utilities/math_ops.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void init_strain_field(kynema_sgf::Field& fld, amrex::Real srate)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    amrex::Real offset =
        (fld.field_location() == kynema_sgf::FieldLoc::CELL) ? 0.5_rt : 0.0_rt;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + offset) * dx[0]);
                const amrex::Real y = problo[1] + ((j + offset) * dx[1]);
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                farrs[nbx](i, j, k, 0) = x / std::sqrt(6.0_rt) * srate;
                farrs[nbx](i, j, k, 1) = y / std::sqrt(6.0_rt) * srate;
                farrs[nbx](i, j, k, 2) = z / std::sqrt(6.0_rt) * srate;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_temperature_field(kynema_sgf::Field& fld, amrex::Real tgrad)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset =
        (fld.field_location() == kynema_sgf::FieldLoc::CELL) ? 0.5_rt : 0.0_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                farrs[nbx](i, j, k, 0) = z * tgrad;
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Shear velocity (S z, 0, w) on level 0, including the ghost cells
void init_shear_velocity(
    kynema_sgf::Field& vel, const amrex::Real srate, const amrex::Real wvel)
{
    const auto& geom = vel.repo().mesh().Geom(0);
    const auto problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize()[2];
    const auto& varrs = vel(0).arrays();
    amrex::ParallelFor(
        vel(0), vel.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            varrs[nbx](i, j, k, 0) = srate * (problo[2] + ((k + 0.5_rt) * dz));
            varrs[nbx](i, j, k, 1) = 0.0_rt;
            varrs[nbx](i, j, k, 2) = wvel;
        });
    amrex::Gpu::streamSynchronize();
}

//! Solid rotation about the vertical axis through (xc, yc) on level 0,
//! including the ghost cells
void init_solid_rotation(
    kynema_sgf::Field& vel,
    const amrex::Real omega,
    const amrex::Real xc,
    const amrex::Real yc)
{
    const auto& geom = vel.repo().mesh().Geom(0);
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const auto& varrs = vel(0).arrays();
    amrex::ParallelFor(
        vel(0), vel.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            varrs[nbx](i, j, k, 0) = -omega * (y - yc);
            varrs[nbx](i, j, k, 1) = omega * (x - xc);
            varrs[nbx](i, j, k, 2) = 0.0_rt;
        });
    amrex::Gpu::streamSynchronize();
}

/** Smallest and largest ratio of the velocity sensor to its exact value
 *
 *  For the velocity (S z, 0, w) the sensor is
 *  -w S^2 z L / (|u| (k + c_u |u|^2)) with the neutral length scale. The
 *  one-sided gradient at the lower and upper boundary is not exact for this
 *  field, so only the interior cells in z are compared.
 *
 *  \param sensor Pressure-gradient sensor field
 *  \param srate Shear rate S
 *  \param wvel Vertical velocity w
 *  \param tke_val Uniform turbulent kinetic energy
 *  \param c_u Velocity weight of the sensor scale
 */
amrex::Array<amrex::Real, 2> velocity_sensor_ratio(
    const kynema_sgf::Field& sensor,
    const amrex::Real srate,
    const amrex::Real wvel,
    const amrex::Real tke_val,
    const amrex::Real c_u)
{
    const auto& geom = sensor.repo().mesh().Geom(0);
    const auto problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize()[2];
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const int klo = geom.Domain().smallEnd(2) + 1;
    const int khi = geom.Domain().bigEnd(2) - 1;
    const amrex::Real huge = std::numeric_limits<amrex::Real>::max();
    const auto& sensor_arrs = sensor(0).const_arrays();
    auto ratio = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMin, amrex::ReduceOpMax>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, sensor(0),
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Real, amrex::Real> {
            if (k < klo || k > khi) {
                return amrex::makeTuple(huge, -huge);
            }
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dz);
            const amrex::Real ux = srate * z;
            const amrex::Real usqr = (ux * ux) + (wvel * wvel);
            const amrex::Real lscale =
                (lambda * kappa * z) / (lambda + (kappa * z));
            const amrex::Real expected =
                -wvel * srate * srate * z * lscale /
                (std::sqrt(usqr) * (tke_val + (c_u * usqr)));
            const amrex::Real r = sensor_arrs[nbx](i, j, k) / expected;
            return amrex::makeTuple(r, r);
        });
    amrex::ParallelDescriptor::ReduceRealMin(amrex::get<0>(ratio));
    amrex::ParallelDescriptor::ReduceRealMax(amrex::get<1>(ratio));
    return {amrex::get<0>(ratio), amrex::get<1>(ratio)};
}

/** Smallest and largest ratio of the eddy viscosity to its value in solid
 *  rotation with the richardson curvature model
 *
 *  With N2_c = 4 Omega^2 the curvature Rt is Rt_c = 4 Omega^2 L^2 / (Cmu^6 k)
 *  and the eddy viscosity is rho Cmu(Rt_c) L sqrt(k), with the neutral length
 *  scale.
 *
 *  \param muturb Eddy viscosity field
 *  \param rho0 Uniform density
 *  \param tke_val Uniform turbulent kinetic energy
 *  \param omega Rotation rate Omega
 */
amrex::Array<amrex::Real, 2> rotation_viscosity_ratio(
    const kynema_sgf::Field& muturb,
    const amrex::Real rho0,
    const amrex::Real tke_val,
    const amrex::Real omega)
{
    const auto& geom = muturb.repo().mesh().Geom(0);
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real Cmu6 = kynema_sgf::utils::powi(Cmu, 6);
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const auto& mu_arrs = muturb(0).const_arrays();
    auto ratio = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMin, amrex::ReduceOpMax>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, muturb(0),
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept
            -> amrex::GpuTuple<amrex::Real, amrex::Real> {
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real lscale =
                (lambda * kappa * z) / (lambda + (kappa * z));
            const amrex::Real rt =
                4.0_rt * omega * omega * lscale * lscale / (Cmu6 * tke_val);
            const amrex::Real cmu_rt =
                (Cmu + (0.108_rt * rt)) /
                (1.0_rt + (0.308_rt * rt) + (0.00837_rt * rt * rt));
            const amrex::Real expected =
                rho0 * cmu_rt * lscale * std::sqrt(tke_val);
            const amrex::Real r = mu_arrs[nbx](i, j, k) / expected;
            return amrex::makeTuple(r, r);
        });
    amrex::ParallelDescriptor::ReduceRealMin(amrex::get<0>(ratio));
    amrex::ParallelDescriptor::ReduceRealMax(amrex::get<1>(ratio));
    return {amrex::get<0>(ratio), amrex::get<1>(ratio)};
}

} // namespace

class TurbRANSTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{10, 10, 64}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{
                {1024.0_rt, 1024.0_rt, 1024.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }

    /** Set up a neutral ABL and create a KLAxell-type turbulence model
     *
     *  \param model Name of the turbulence model
     */
    void create_klaxell_model(const std::string& model);

    /** Check the eddy viscosity of a KLAxell-type model on a uniform strain
     *  field
     *
     *  \param model Name of the turbulence model
     */
    void check_klaxell_viscosity(const std::string& model);

    //! Density of the test case
    const amrex::Real m_rho0{1.2_rt};

    //! Reference temperature of the test case
    const amrex::Real m_tref{265.0_rt};
};

void TurbRANSTest::create_klaxell_model(const std::string& model)
{
    // Parser inputs for turbulence model
    const amrex::Real gravz = 10.0_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", model);
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", m_rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -gravz};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        pp.add("initial_wind_profile", true);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        pp.addarr("wind_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{m_tref, m_tref, m_tref};
        pp.addarr("temperature_values", t_vals);
        amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
        pp.addarr("u_values", u_vals);
        amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("v_values", v_vals);
        amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
        pp.addarr("tke_values", tke_vals);
        pp.add("surface_temp_flux", 0.0_rt);
    }
    // Transport
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", m_tref);
    }

    // Initialize necessary parts of solver
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Create turbulence model
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
}

void TurbRANSTest::check_klaxell_viscosity(const std::string& model)
{
    create_klaxell_model(model);
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    // Constants for fields
    const amrex::Real srate = 20.0_rt;
    const amrex::Real Tgz = 0.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real x3 = 1016.0_rt;
    const amrex::Real lscale_s =
        (lambda * kappa * x3) / (lambda + (kappa * x3));
    const amrex::Real tlscale_val = lscale_s;
    const amrex::Real tke_val = 0.1_rt;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_strain_field(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(m_rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    init_temperature_field(temp, Tgz);
    // Give values to tlscale and tke arrays
    auto& tlscale = sim().repo().get_field("turb_lscale");
    tlscale.setVal(tlscale_val);
    auto& tke = sim().repo().get_field("tke");
    tke.setVal(tke_val);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto max_val = utils::field_max(muturb);
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real epsilon = kynema_sgf::utils::powi(Cmu, 3) *
                                std::pow(tke_val, 1.5_rt) /
                                (tlscale_val + 1.0e-3_rt);
    const amrex::Real stratification = 0.0_rt;
    const amrex::Real Rt =
        kynema_sgf::utils::powi(tke_val / epsilon, 2) * stratification;
    const amrex::Real Cmu_Rt = (0.556_rt + (0.108_rt * Rt)) /
                               (1.0_rt + (0.308_rt * Rt) +
                                (0.00837_rt * kynema_sgf::utils::powi(Rt, 2)));
    const amrex::Real tol = 0.12_rt;
    const amrex::Real nut_max =
        m_rho0 * Cmu_Rt * tlscale_val * std::sqrt(tke_val);
    EXPECT_NEAR(max_val, nut_max, tol);
}

TEST_F(TurbRANSTest, test_1eqKrans_setup_calc)
{
    check_klaxell_viscosity("KLAxell");
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_setup_calc)
{
    check_klaxell_viscosity("KLAxellSeparation");
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_gate_default)
{
    // The gate relaxes over 10 s by default, but without the sensor there is
    // no gate and the default does not require the sensor
    create_klaxell_model("KLAxellSeparation");
    EXPECT_EQ(
        sim().turbulence_model().model_coeffs().at("gate_relaxation_time"),
        10.0_rt);
    EXPECT_FALSE(sim().repo().field_exists("separation_gate"));
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_gate_default_sensor)
{
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("pressure_gradient_sensor", true);
    }
    create_klaxell_model("KLAxellSeparation");
    EXPECT_TRUE(sim().repo().field_exists("separation_gate"));
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_gate_requires_sensor)
{
    // A relaxation time set in the input without the sensor is an error
    {
        amrex::ParmParse pp("KLAxellSeparation_coeffs");
        pp.add("gate_relaxation_time", 10.0_rt);
    }
    EXPECT_THROW(create_klaxell_model("KLAxellSeparation"), std::runtime_error);
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_negative_strength)
{
    // A negative limiter strength could make the limiter denominator vanish
    {
        amrex::ParmParse pp("KLAxellSeparation_coeffs");
        pp.add("realizable_cmu_strength", -1.0_rt);
    }
    EXPECT_THROW(create_klaxell_model("KLAxellSeparation"), std::runtime_error);
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_negative_velocity_weight)
{
    // A negative velocity weight could make the sensor scale vanish
    {
        amrex::ParmParse pp("KLAxellSeparation_coeffs");
        pp.add("sensor_velocity_weight", -0.05_rt);
    }
    EXPECT_THROW(create_klaxell_model("KLAxellSeparation"), std::runtime_error);
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_pressure_gradient_sensor)
{
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("pressure_gradient_sensor", true);
    }
    create_klaxell_model("KLAxellSeparation");
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();

    const amrex::Real tke_val = 0.1_rt;
    const amrex::Real gp_val = 0.3_rt;
    auto& vel = repo.get_field("velocity");
    auto& temp = repo.get_field("temperature");
    auto& gp = repo.get_field("gp");
    const auto& sensor = repo.get_field("pressure_gradient_sensor");
    repo.get_field("density").setVal(m_rho0);
    repo.get_field("tke").setVal(tke_val);
    repo.get_field("turb_lscale").setVal(1.0_rt);
    temp.setVal(amrex::Vector<amrex::Real>{m_tref}, temp.num_grow()[0]);

    const auto update = [&](const amrex::Vector<amrex::Real>& uvec,
                            const amrex::Vector<amrex::Real>& gpvec) {
        vel.setVal(uvec, vel.num_grow()[0]);
        gp.setVal(gpvec);
        tmodel.update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    };

    // The sensor is (u / |u|) . grad(p) L / (rho (k + c_u |u|^2)). In a
    // neutral flow the length scale is lambda kappa z / (lambda + kappa z);
    // the lowest and highest cell centers are at 8 m and 1016 m.
    const amrex::Real c_u = tmodel.model_coeffs().at("sensor_velocity_weight");
    EXPECT_EQ(c_u, 0.05_rt);
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const auto lscale = [=](const amrex::Real z) {
        return (lambda * kappa * z) / (lambda + (kappa * z));
    };
    const amrex::Real lscale_bottom = lscale(8.0_rt);
    const amrex::Real lscale_top = lscale(1016.0_rt);
    const amrex::Real rtol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    // Check the sensor range for a sensor per unit length scale
    const auto check_range = [&](const amrex::Real per_length) {
        EXPECT_NEAR(
            utils::field_min(sensor), per_length * lscale_bottom,
            rtol * per_length * lscale_bottom);
        EXPECT_NEAR(
            utils::field_max(sensor), per_length * lscale_top,
            rtol * per_length * lscale_top);
    };

    // Pressure gradient along the flow
    update({8.0_rt, 0.0_rt, 0.0_rt}, {gp_val, 0.0_rt, 0.0_rt});
    check_range(gp_val / (m_rho0 * (tke_val + (c_u * 64.0_rt))));

    // Pressure gradient normal to the flow: no signal
    update({8.0_rt, 0.0_rt, 0.0_rt}, {0.0_rt, gp_val, 0.0_rt});
    EXPECT_EQ(utils::field_min(sensor), 0.0_rt);
    EXPECT_EQ(utils::field_max(sensor), 0.0_rt);

    // Vertical pressure gradient with a vertical velocity: the sensor keeps
    // the vertical term
    update({8.0_rt, 0.0_rt, 4.0_rt}, {0.0_rt, 0.0_rt, gp_val});
    check_range(
        4.0_rt / std::sqrt(80.0_rt) * gp_val /
        (m_rho0 * (tke_val + (c_u * 80.0_rt))));

    // Still air: the sensor stays finite and is zero
    update({0.0_rt, 0.0_rt, 0.0_rt}, {gp_val, 0.0_rt, 0.0_rt});
    EXPECT_EQ(utils::field_min(sensor), 0.0_rt);
    EXPECT_EQ(utils::field_max(sensor), 0.0_rt);

    // Nearly vanishing TKE: the velocity term keeps the sensor bounded
    const amrex::Real tke_small = 1.0e-6_rt;
    repo.get_field("tke").setVal(tke_small);
    update({8.0_rt, 0.0_rt, 0.0_rt}, {gp_val, 0.0_rt, 0.0_rt});
    check_range(gp_val / (m_rho0 * (tke_small + (c_u * 64.0_rt))));
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_velocity_sensor)
{
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("pressure_gradient_sensor", true);
        pp.add("sensor_source", std::string("velocity"));
    }
    create_klaxell_model("KLAxellSeparation");
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();

    const amrex::Real tke_val = 0.1_rt;
    const amrex::Real srate = 0.01_rt;
    auto& vel = repo.get_field("velocity");
    auto& temp = repo.get_field("temperature");
    const auto& sensor = repo.get_field("pressure_gradient_sensor");
    repo.get_field("density").setVal(m_rho0);
    repo.get_field("tke").setVal(tke_val);
    repo.get_field("turb_lscale").setVal(1.0_rt);
    temp.setVal(amrex::Vector<amrex::Real>{m_tref}, temp.num_grow()[0]);
    // A pressure gradient that the velocity sensor must ignore
    repo.get_field("gp").setVal(
        amrex::Vector<amrex::Real>{1.0e3_rt, 1.0e3_rt, 1.0e3_rt});

    // Velocity (S z, 0, w), including the ghost cells. The device kernels
    // live in free functions: nvcc rejects them in a test body.
    const auto update = [&](const amrex::Real wvel) {
        init_shear_velocity(vel, srate, wvel);
        tmodel.update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    };

    // Shear without change along the flow: no signal
    update(0.0_rt);
    EXPECT_EQ(utils::field_min(sensor), 0.0_rt);
    EXPECT_EQ(utils::field_max(sensor), 0.0_rt);

    // Flow down into slower air: (u . grad)(|u|^2 / 2) = w S^2 z, so the
    // sensor is -w S^2 z L / (|u| (k + c_u |u|^2)) with the neutral length
    // scale. The one-sided gradient at the lower and upper boundary is not
    // exact for this field, so only the interior cells are checked.
    const amrex::Real wvel = -2.0_rt;
    update(wvel);
    const amrex::Real c_u = tmodel.model_coeffs().at("sensor_velocity_weight");
    const auto ratio = velocity_sensor_ratio(sensor, srate, wvel, tke_val, c_u);
    const amrex::Real rtol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    EXPECT_NEAR(ratio[0], 1.0_rt, rtol);
    EXPECT_NEAR(ratio[1], 1.0_rt, rtol);
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_realizable_cmu)
{
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("pressure_gradient_sensor", true);
        pp.add("realizable_cmu", true);
    }
    {
        // Instantaneous gate, so that one update applies the limiter
        amrex::ParmParse pp("KLAxellSeparation_coeffs");
        pp.add("gate_relaxation_time", 0.0_rt);
    }
    create_klaxell_model("KLAxellSeparation");
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();
    const auto coeffs = tmodel.model_coeffs();
    EXPECT_EQ(coeffs.at("sensor_threshold"), 0.1_rt);
    EXPECT_EQ(coeffs.at("realizable_cmu_strength"), 1.0_rt);

    // Uniform strain rate S, uniform TKE and a stable temperature gradient.
    // Below meso_sponge_start with no surface flux the neutral length scale
    // applies, so Rt = 0 and Cmu(Rt) = Cmu.
    const amrex::Real srate = 1.0_rt;
    const amrex::Real tke_val = 0.1_rt;
    const amrex::Real tgrad = 0.01_rt;
    init_strain_field(repo.get_field("velocity"), srate);
    init_temperature_field(repo.get_field("temperature"), tgrad);
    repo.get_field("density").setVal(m_rho0);
    repo.get_field("tke").setVal(tke_val);
    repo.get_field("turb_lscale").setVal(1.0_rt);
    auto& gp = repo.get_field("gp");
    const auto& sensor = repo.get_field("pressure_gradient_sensor");
    const auto& muturb = repo.get_field("mu_turb");
    const auto& shear_prod = repo.get_field("shear_prod");
    const auto& buoy_prod = repo.get_field("buoy_prod");
    const auto update = [&](const amrex::Real gp_val) {
        gp.setVal(amrex::Vector<amrex::Real>{gp_val, gp_val, gp_val});
        tmodel.update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    };
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real rtol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // No pressure gradient: the sensor stays below the threshold and the eddy
    // viscosity keeps the KLAxell value rho Cmu L sqrt(k), largest at the top
    update(0.0_rt);
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real z_top = 1016.0_rt;
    const amrex::Real lscale_top =
        (lambda * kappa * z_top) / (lambda + (kappa * z_top));
    const amrex::Real mu_top = m_rho0 * Cmu * lscale_top * std::sqrt(tke_val);
    EXPECT_NEAR(utils::field_max(muturb), mu_top, rtol * mu_top);
    // The buoyancy production is -mu N^2 with the same eddy viscosity
    const amrex::Real nsqr =
        -utils::field_min(buoy_prod) / utils::field_max(muturb);
    EXPECT_GT(nsqr, 0.0_rt);

    // Strong pressure gradient along the flow: the sensor fires everywhere,
    // Sigma = L S / sqrt(k) exceeds Cmu in every cell, and with c_s = 1 the
    // eddy viscosity becomes rho Cmu^2 k / S
    update(1.0e4_rt);
    EXPECT_GT(utils::field_min(sensor), coeffs.at("sensor_threshold"));
    const amrex::Real mu_limited = m_rho0 * Cmu * Cmu * tke_val / srate;
    EXPECT_NEAR(utils::field_min(muturb), mu_limited, rtol * mu_limited);
    EXPECT_NEAR(utils::field_max(muturb), mu_limited, rtol * mu_limited);
    // The shear and buoyancy production follow the limited eddy viscosity
    const amrex::Real shear_limited = srate * srate * mu_limited;
    EXPECT_NEAR(
        utils::field_min(shear_prod), shear_limited, rtol * shear_limited);
    EXPECT_NEAR(
        utils::field_max(shear_prod), shear_limited, rtol * shear_limited);
    const amrex::Real buoy_limited = -nsqr * mu_limited;
    EXPECT_NEAR(
        utils::field_min(buoy_prod), buoy_limited,
        rtol * std::abs(buoy_limited));
    EXPECT_NEAR(
        utils::field_max(buoy_prod), buoy_limited,
        rtol * std::abs(buoy_limited));
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_gate_relaxation)
{
    const amrex::Real tau = 10.0_rt;
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("pressure_gradient_sensor", true);
        pp.add("realizable_cmu", true);
    }
    {
        amrex::ParmParse pp("KLAxellSeparation_coeffs");
        pp.add("gate_relaxation_time", tau);
    }
    create_klaxell_model("KLAxellSeparation");
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();
    EXPECT_EQ(tmodel.model_coeffs().at("gate_relaxation_time"), tau);

    // Uniform strain rate and TKE with a neutral temperature field, so
    // Cmu(Rt) = Cmu and Sigma = L S / sqrt(k) exceeds Cmu in every cell
    const amrex::Real srate = 1.0_rt;
    const amrex::Real tke_val = 0.1_rt;
    init_strain_field(repo.get_field("velocity"), srate);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("density").setVal(m_rho0);
    repo.get_field("tke").setVal(tke_val);
    repo.get_field("turb_lscale").setVal(1.0_rt);
    auto& gp = repo.get_field("gp");
    const auto& sensor = repo.get_field("pressure_gradient_sensor");
    const auto& muturb = repo.get_field("mu_turb");
    const auto& gate = repo.get_field("separation_gate");
    const amrex::Real dt = 1.0_rt;
    sim().time().delta_t() = dt;
    // One time step: viscosity update with the stored gate, then relaxation
    const auto step = [&](const amrex::Real gp_val) {
        gp.setVal(amrex::Vector<amrex::Real>{gp_val, gp_val, gp_val});
        tmodel.update_turbulent_viscosity(
            kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
        tmodel.post_advance_work();
    };
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const auto lscale = [=](const amrex::Real z) {
        return (lambda * kappa * z) / (lambda + (kappa * z));
    };
    // Eddy viscosity for a gate g; it grows with L, so the extremes are in
    // the lowest and highest cells (8 m and 1016 m)
    const auto mu_gated = [=, this](const amrex::Real g, const amrex::Real z) {
        const amrex::Real sigma = lscale(z) * srate / std::sqrt(tke_val);
        return m_rho0 * Cmu * lscale(z) * std::sqrt(tke_val) /
               (1.0_rt + (g * ((sigma / Cmu) - 1.0_rt)));
    };
    const amrex::Real rtol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    const auto check_gate = [&](const amrex::Real g) {
        EXPECT_NEAR(utils::field_min(gate), g, rtol);
        EXPECT_NEAR(utils::field_max(gate), g, rtol);
    };
    EXPECT_EQ(utils::field_max(gate), 0.0_rt);

    // First step with a sensor far above twice the threshold: the viscosity
    // still uses the initial gate 0, and the gate relaxes toward 1
    const amrex::Real gp_strong = 1.0e5_rt;
    step(gp_strong);
    EXPECT_GT(
        utils::field_min(sensor),
        2.0_rt * tmodel.model_coeffs().at("sensor_threshold"));
    EXPECT_NEAR(
        utils::field_max(muturb), mu_gated(0.0_rt, 1016.0_rt),
        rtol * mu_gated(0.0_rt, 1016.0_rt));
    const amrex::Real decay = std::exp(-dt / tau);
    const amrex::Real g1 = 1.0_rt - decay;
    check_gate(g1);

    // Second step: the viscosity uses g1
    step(gp_strong);
    EXPECT_NEAR(
        utils::field_min(muturb), mu_gated(g1, 8.0_rt),
        rtol * mu_gated(g1, 8.0_rt));
    EXPECT_NEAR(
        utils::field_max(muturb), mu_gated(g1, 1016.0_rt),
        rtol * mu_gated(g1, 1016.0_rt));
    const amrex::Real g2 = 1.0_rt - (decay * decay);
    check_gate(g2);

    // No pressure gradient: the gate decays toward 0
    step(0.0_rt);
    check_gate(g2 * decay);
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_curvature_richardson)
{
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("curvature_correction", true);
        pp.add("curvature_model", std::string("richardson"));
    }
    create_klaxell_model("KLAxellSeparation");
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();
    EXPECT_EQ(tmodel.model_coeffs().at("curvature_coefficient"), 1.0_rt);

    const amrex::Real rho0 = m_rho0;
    const amrex::Real tke_val = 0.1_rt;
    auto& vel = repo.get_field("velocity");
    auto& temp = repo.get_field("temperature");
    const auto& muturb = repo.get_field("mu_turb");
    repo.get_field("density").setVal(rho0);
    repo.get_field("tke").setVal(tke_val);
    repo.get_field("turb_lscale").setVal(1.0_rt);
    temp.setVal(amrex::Vector<amrex::Real>{m_tref}, temp.num_grow()[0]);

    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real rtol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Simple shear u = (S z, 0, 0): straight streamlines keep the KLAxell
    // eddy viscosity rho Cmu L sqrt(k), largest at the top
    init_shear_velocity(vel, 0.01_rt, 0.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const amrex::Real z_top = 1016.0_rt;
    const amrex::Real lscale_top =
        (lambda * kappa * z_top) / (lambda + (kappa * z_top));
    const amrex::Real mu_top = rho0 * Cmu * lscale_top * std::sqrt(tke_val);
    EXPECT_NEAR(utils::field_max(muturb), mu_top, rtol * mu_top);

    // Solid rotation about a vertical axis: N2_c = 4 Omega^2, so
    // Rt_c = 4 Omega^2 L^2 / (Cmu^6 k) and mu = rho Cmu(Rt_c) L sqrt(k). The
    // axis lies between cell centers, so |u| > 0 everywhere.
    const amrex::Real omega = 1.0e-3_rt;
    init_solid_rotation(vel, omega, 512.0_rt, 512.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto ratio = rotation_viscosity_ratio(muturb, rho0, tke_val, omega);
    EXPECT_NEAR(ratio[0], 1.0_rt, rtol);
    EXPECT_NEAR(ratio[1], 1.0_rt, rtol);
}

TEST_F(TurbRANSTest, test_1eqKrans_separation_curvature_rotation_function)
{
    {
        amrex::ParmParse pp("KLAxellSeparation");
        pp.add("curvature_correction", true);
        pp.add("curvature_model", std::string("rotation_function"));
    }
    create_klaxell_model("KLAxellSeparation");
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();

    const amrex::Real tke_val = 0.1_rt;
    auto& vel = repo.get_field("velocity");
    auto& temp = repo.get_field("temperature");
    const auto& muturb = repo.get_field("mu_turb");
    repo.get_field("density").setVal(m_rho0);
    repo.get_field("tke").setVal(tke_val);
    repo.get_field("turb_lscale").setVal(1.0_rt);
    temp.setVal(amrex::Vector<amrex::Real>{m_tref}, temp.num_grow()[0]);

    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real rtol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Simple shear: S = Omega, so r* = 1, r~ = 0 and f = 1; the eddy
    // viscosity keeps the KLAxell value
    init_shear_velocity(vel, 0.01_rt, 0.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const amrex::Real z_top = 1016.0_rt;
    const amrex::Real lscale_top =
        (lambda * kappa * z_top) / (lambda + (kappa * z_top));
    const amrex::Real mu_top = m_rho0 * Cmu * lscale_top * std::sqrt(tke_val);
    EXPECT_NEAR(utils::field_max(muturb), mu_top, rtol * mu_top);

    // Solid rotation: S = 0, so r* = 0 and f = -c_r1, limited to 0; the eddy
    // viscosity vanishes
    init_solid_rotation(vel, 1.0e-3_rt, 512.0_rt, 512.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    EXPECT_EQ(utils::field_min(muturb), 0.0_rt);
    EXPECT_EQ(utils::field_max(muturb), 0.0_rt);
}

} // namespace kynema_sgf_tests
