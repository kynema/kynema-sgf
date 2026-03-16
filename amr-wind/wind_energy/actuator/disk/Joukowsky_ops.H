#ifndef JOUKOWSKY_OPS_H
#define JOUKOWSKY_OPS_H

#include <algorithm>

#include "amr-wind/wind_energy/actuator/disk/Joukowsky.H"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops {

namespace joukowsky {
void parse_and_gather_params(const utils::ActParser& pp, JoukowskyData& data);
AMREX_INLINE void
set_current_angular_velocity(JoukowskyData& data, const amrex::Real uInfSqr)
{
    const amrex::Real uInfMag = std::sqrt(uInfSqr);
    data.current_angular_velocity = ::amr_wind::interp::linear(
        data.table_velocity, data.angular_velocity, uInfMag);
}
void prepare_netcdf_file(
    const std::string& name,
    const JoukowskyData& data,
    const ActInfo& info,
    const ActGrid& grid);

void write_netcdf(
    const std::string& name,
    const JoukowskyData& data,
    const ActInfo& info,
    const ActGrid& /*unused*/,
    amrex::Real time);

void update_disk_points(Joukowsky::DataType& data);

} // namespace joukowsky

template <>
struct ReadInputsOp<Joukowsky, ActSrcDisk>
{
    void operator()(Joukowsky::DataType& data, const utils::ActParser& pp)
    {
        auto& meta = data.meta();
        joukowsky::parse_and_gather_params(pp, meta);
        base::do_parse_based_computations<Joukowsky>(data);
    }
};

template <>
struct UpdatePosOp<Joukowsky, ActSrcDisk>
{
    void operator()(typename Joukowsky::DataType& data)
    {

        if (!data.sim().helics().is_activated()) {
            return;
        }

#ifdef AMR_WIND_USE_HELICS

        auto& meta = data.meta();

        const amrex::Real turbine_direction =
            -data.sim().helics().m_turbine_yaw_to_amrwind[data.info().id] +
            90.0_rt;
        const amrex::Real turbine_direction_radian =
            amr_wind::utils::radians(turbine_direction);

        meta.normal_vec[0] = std::cos(turbine_direction_radian);
        meta.normal_vec[1] = std::sin(turbine_direction_radian);
        meta.normal_vec[2] = 0.0_rt;

        meta.sample_vec[0] = meta.normal_vec[0];
        meta.sample_vec[1] = meta.normal_vec[1];
        meta.sample_vec[2] = meta.normal_vec[2];

        joukowsky::update_disk_points(data);
#endif
    }
};

template <>
struct InitDataOp<Joukowsky, ActSrcDisk>
{
    void operator()(typename Joukowsky::DataType& data)
    {
        ops::base::allocate_basic_grid_quantities<Joukowsky>(data);

        auto& meta = data.meta();
        meta.tip_correction.resize(meta.num_vel_pts_r);
        meta.root_correction.resize(meta.num_vel_pts_r);
        meta.f_normal.resize(meta.num_vel_pts_r);
        meta.f_theta.resize(meta.num_force_pts);

        std::ranges::fill(meta.tip_correction, 1.0_rt);

        if (meta.use_root_correction) {
            const amrex::Real dx = 1.0_rt / meta.num_vel_pts_r;

            const auto a = meta.root_correction_coefficient;
            const auto b = meta.root_correction_exponent;

            const auto delta = amrex::max<amrex::Real>(
                meta.vortex_core_size / meta.radius(),
                vs::DTraits<amrex::Real>::eps());

            const auto factor = dx / delta;

            for (int i = 0; i < meta.root_correction.size(); ++i) {
                meta.root_correction[i] =
                    1.0_rt - std::exp(-a * std::pow((i + 0.5_rt) * factor, b));
            }
        } else {
            std::ranges::fill(meta.root_correction, 1.0_rt);
        }

        joukowsky::update_disk_points(data);
    }
};

/**
  \brief Compute Forces following procedures of Sorenson 2020

  Follow the procedure for computing body forces and applying tip/root
  corrections using a Joukowsky disk as outlined in:

  Sørensen, Jens Nørkær, et al. "Analytical body forces in numerical actuator
  disc model of wind turbines." Renewable Energy 147 (2020): 2259-2271.

  Upgraded to handle rated turbine conditions following the implementation
  Sørensen, J. N.: Generalized analytical body force model for actuator
  disc computations of wind turbines,
  Wind Energ. Sci. Discuss. [preprint], https://doi.org/10.5194/wes-2022-108
  https://wes.copernicus.org/preprints/wes-2022-108/
 */
template <>
struct ComputeForceOp<Joukowsky, ActSrcDisk>
{
    void operator()(Joukowsky::DataType& data)
    {
        // Equation comments refer to Sorenson 2020 (full reference above)
        auto& grid = data.grid();
        auto& ddata = data.meta();
        const amrex::Real machine_eps = vs::DTraits<amrex::Real>::eps();
        ddata.disk_force *= 0.0_rt;

        const amrex::Real uInfSqr = base::compute_reference_velocity_sqr(ddata);
        const amrex::Real U_ref =
            amrex::max<amrex::Real>(std::sqrt(uInfSqr), machine_eps);
        base::set_thrust_coefficient(ddata, uInfSqr);
        joukowsky::set_current_angular_velocity(ddata, uInfSqr);

        const amrex::Real Ct = ddata.current_ct;
        const amrex::Real dx = ddata.dr / ddata.radius();

        const ops::base::AreaComputer area(
            ddata.radius(), ddata.num_vel_pts_r, ddata.num_vel_pts_t);

        const amrex::Real lambda =
            ddata.radius() * ddata.current_angular_velocity / U_ref;

        ddata.current_tip_speed_ratio = lambda;

        const amrex::Real lambda_2 = lambda * lambda;

        // step 0: compute tip correction based on current wind speed (eq 9)
        if (ddata.use_tip_correction) {
            auto& f = ddata.tip_correction;
            for (int i = 0; i < f.size(); ++i) {
                const auto x = (i + 0.5_rt) * dx;
                f[i] = 2.0_rt / std::numbers::pi_v<amrex::Real> *
                       std::acos(
                           std::exp(
                               -0.5_rt * ddata.num_blades *
                               std::sqrt(1.0_rt + lambda_2) * (1.0_rt - x)));
            }
        }

        // Compute the solid body term S0 (Sorensen 2022, eq. 27)
        amrex::Real S0 = 0.0_rt;
        const amrex::Real alpha1 = ddata.S0_alpha1;
        const amrex::Real alpha2 = ddata.S0_alpha2;

        if (ddata.Ct_rated > 0.0_rt) {
            if (Ct < ddata.Ct_rated) {
                S0 = alpha1 *
                     std::pow((ddata.Ct_rated - Ct) / ddata.Ct_rated, 3.0_rt);
            } else {
                S0 = alpha2 * (ddata.Ct_rated - Ct) / ddata.Ct_rated;
            }
        }

        // step 1: compute Ct, a1 and a2 coefficients (eq 16)
        // a3, a4, a5 coefficients are from Sorensen, 2022
        amrex::Real a1 = 0.0_rt;
        amrex::Real a2 = 0.0_rt;
        amrex::Real a3 = 0.0_rt;
        amrex::Real a4 = 0.0_rt;
        amrex::Real a5 = 0.0_rt;

        for (int ip = 0; ip < ddata.num_vel_pts_r; ip++) {
            const amrex::Real x =
                amrex::max<amrex::Real>((ip + 0.5_rt) * dx, machine_eps);

            a1 += std::pow(ddata.tip_correction[ip], 2.0_rt) *
                  std::pow(ddata.root_correction[ip], 2.0_rt) / x * dx;

            a2 += ddata.tip_correction[ip] * ddata.root_correction[ip] * x * dx;

            a3 += std::pow(ddata.tip_correction[ip], 2.0_rt) *
                  std::pow(ddata.root_correction[ip], 2.0_rt) * x * dx;

            a4 += ddata.tip_correction[ip] * ddata.root_correction[ip] *
                  std::pow(x, 3.0_rt) * dx;

            a5 += std::pow(ddata.tip_correction[ip], 2.0_rt) *
                  std::pow(ddata.root_correction[ip], 2.0_rt) *
                  std::pow(x, 3.0_rt) * dx;
        }

        // step 2: determine the circulation (q0) from a1, a2 and Ct
        //         (eq. 18, Sorensen, 2022)

        //      sqrt[(a2*L+a3*S0)^2+a1*(0.5_rt*Ct-2*a4*L*S0-a5*S0^2)]-(a2*L+a3*S0)
        // q0 = --------------------------------------------------------------
        //                            a1

        const amrex::Real term1 = std::pow((a2 * lambda) + (a3 * S0), 2.0_rt);
        const amrex::Real term2 =
            a1 * (0.5_rt * Ct - 2.0_rt * a4 * lambda * S0 - a5 * S0 * S0);
        const amrex::Real term3 = (a2 * lambda) + (a3 * S0);

        const amrex::Real numerator = std::sqrt(term1 + term2) - term3;

        const amrex::Real denominator =
            amrex::max<amrex::Real>(a1, machine_eps);

        const amrex::Real q0 = numerator / denominator;

        // step 3: compute normal force (fz) and azimuthal force (f_theta) (eq
        // 13) and cp (eq 20)

        VecSlice disk_velocity = ::amr_wind::utils::slice(
            grid.vel, ddata.num_force_pts, ddata.num_force_pts);

        amrex::Real moment = 0.0_rt;

        int ip = 0;

        const auto dTheta = amr_wind::utils::two_pi() / ddata.num_vel_pts_t;
        const auto current_time = data.sim().time().current_time();

        for (int i = 0; i < ddata.num_vel_pts_r; i++) {
            const amrex::Real& F = ddata.tip_correction[i];
            const amrex::Real& g = ddata.root_correction[i];

            const amrex::Real x =
                amrex::max<amrex::Real>((i + 0.5_rt) * dx, machine_eps);

            const amrex::Real f_z =
                (q0 / x + S0 * x) * g * F *
                (lambda * x + 0.5_rt * (q0 / x + S0 * x) * g * F);

            ddata.f_normal[i] = f_z;

            const auto dArea = area.area_section(i);

            for (int j = 0; j < ddata.num_vel_pts_t; j++, ip++) {
                const auto point_current = grid.pos[ip];
                const auto theta_angle = j * dTheta;

                // TODO add sign convention for rotation (+/- RPM)
                // Negative sign on theta_vec means turbine is rotating
                // clockwise when standing upstream looking downstream
                const auto theta_vec = -utils::compute_tangential_vector(
                    ddata.center, ddata.normal_vec, point_current);

                // normal vec by definition is opposite of the wind direction
                // so we need to flip the sign to give the actual disk velocity
                const amrex::Real u_disk_ij =
                    -ddata.normal_vec & disk_velocity[ip];

                const amrex::Real f_theta =
                    u_disk_ij / U_ref * (q0 / x + S0 * x) * g * F;
                ddata.f_theta[ip] = f_theta;

                // eq 18
                moment += x * ddata.radius() * f_theta * dArea;

                const amrex::Real f_awc =
                    ddata.awc_amplitude *
                    std::cos(
                        (ddata.awc_angular_frequency * current_time) -
                        (ddata.awc_azimuthal_mode * theta_angle) +
                        ddata.awc_clocking_angle);

                grid.force[ip] = ((1 + f_awc) * f_z * ddata.normal_vec +
                                  f_theta * theta_vec) *
                                 ddata.density * uInfSqr * dArea;

                ddata.disk_force = ddata.disk_force + grid.force[ip];
            }
        }

        // equation 20
        const amrex::Real eq_20_denominator = amrex::max<amrex::Real>(
            0.5_rt * std::numbers::pi_v<amrex::Real> * ddata.density *
                std::pow(ddata.radius(), 2.0_rt) * std::pow(U_ref, 3.0_rt),
            machine_eps);

        ddata.current_power =
            ddata.density * uInfSqr * moment * ddata.current_angular_velocity;

        ddata.current_cp = ddata.density * uInfSqr * moment *
                           ddata.current_angular_velocity / eq_20_denominator;
#ifdef AMR_WIND_USE_HELICS
        if (data.info().is_root_proc && data.sim().helics().is_activated()) {

            std::cout << "turbine: " << data.info().id
                      << " jouk power: " << ddata.current_power
                      << " uinfsqr: " << uInfSqr << " moment: " << moment
                      << " ang vel: " << ddata.current_angular_velocity
                      << " normal vec: " << ddata.normal_vec[0] << ' '
                      << ddata.normal_vec[1] << std::endl;

            data.sim().helics().m_turbine_power_to_controller[data.info().id] =
                ddata.current_power;
            const amrex::Real turbine_angle = std::atan2(
                ddata.reference_velocity[1], ddata.reference_velocity[0]);
            data.sim()
                .helics()
                .m_turbine_wind_direction_to_controller[data.info().id] =
                270.0_rt - amr_wind::utils::degrees(turbine_angle);
        }
#endif
    }
};

template <>
struct ProcessOutputsOp<Joukowsky, ActSrcDisk>
{
private:
    // cppcheck-suppress uninitMemberVarPrivate
    const Joukowsky::DataType& m_data;
    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir;

    //! NetCDF output filename for this turbine
    std::string m_nc_filename;

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    // cppcheck-suppress constParameter
    explicit ProcessOutputsOp<Joukowsky, ActSrcDisk>(
        const Joukowsky::DataType& data)
        : m_data(data)
    {}
    void operator()(Joukowsky::DataType& /*data*/) {}
    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }
    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        joukowsky::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }
    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        joukowsky::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace amr_wind::actuator::ops

#endif /* Joukowsky_OPS_H */
