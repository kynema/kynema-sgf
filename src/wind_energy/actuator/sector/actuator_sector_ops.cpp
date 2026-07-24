#include "src/wind_energy/actuator/sector/actuator_sector_ops.H"

#include "src/CFDSim.H"
#include "src/core/Field.H"
#include "src/core/gpu_utils.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/linear_interpolation.H"
#include "src/utilities/ncutils/nc_interface.H"
#include "src/utilities/trig_ops.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "AMReX_Print.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::sector {
namespace {

amrex::Real max_interp_chord(const ActuatorSectorData& meta)
{
    return *std::ranges::max_element(meta.chord_inp);
}

amrex::Real
local_epsilon(const ActuatorSectorData& meta, const amrex::Real chord)
{
    if (meta.user_epsilon) {
        return meta.epsilon;
    }
    return amrex::max(meta.epsilon_min, meta.epsilon_chord * chord);
}

amrex::Real interp_metric_to_radius(
    const amrex::Real target_metric,
    const RealList& metric,
    const RealList& radius)
{
    if (target_metric <= metric.front()) {
        return radius.front();
    }
    if (target_metric >= metric.back()) {
        return radius.back();
    }
    const auto hi =
        std::upper_bound(metric.begin(), metric.end(), target_metric);
    const int ihi = static_cast<int>(std::distance(metric.begin(), hi));
    const int ilo = ihi - 1;
    const amrex::Real t =
        (target_metric - metric[ilo]) / (metric[ihi] - metric[ilo]);
    return (1.0_rt - t) * radius[ilo] + t * radius[ihi];
}

void require_midpoint(const std::string& value, const std::string& name)
{
    if (amrex::toLower(value) != "midpoint") {
        amrex::Abort(
            "ActuatorSector: " + name + " currently supports only 'midpoint'");
    }
}

} // namespace

amrex::Real timestep_width(const CFDSim& sim)
{
    const auto& time = sim.time();
    amrex::Real dt = time.new_time() - time.current_time();
    if (dt <= 0.0_rt) {
        dt = time.delta_t();
    }
    return amrex::max(dt, 0.0_rt);
}

void build_radial_grid(ActuatorSectorData& meta)
{
    const int dense_n = 20000;
    RealList dense_r(dense_n);
    RealList metric(dense_n, 0.0_rt);
    const amrex::Real span = meta.rotor_radius - meta.root_radius;
    AMREX_ALWAYS_ASSERT(span > 0.0_rt);

    for (int i = 0; i < dense_n; ++i) {
        const amrex::Real xi =
            static_cast<amrex::Real>(i) / static_cast<amrex::Real>(dense_n - 1);
        dense_r[i] = meta.root_radius + xi * span;
    }

    // Build a cumulative spacing metric whose local density enforces both
    // requested Gaussian resolution (epsilon_dr / epsilon) and the minimum
    // blade-shape resolution (min_chord_dr / chord). The final station count is
    // the integrated metric, so users specify resolution targets rather than a
    // fixed number of radial stations.
    for (int i = 1; i < dense_n; ++i) {
        const auto density_at = [&](const amrex::Real r) {
            const amrex::Real mu =
                (r - meta.root_radius) / (meta.rotor_radius - meta.root_radius);
            const amrex::Real chord = ::kynema_sgf::interp::linear(
                meta.span_locs, meta.chord_inp, mu);
            const amrex::Real eps = local_epsilon(meta, chord);
            return amrex::max(meta.epsilon_dr / eps, meta.min_chord_dr / chord);
        };
        const amrex::Real d0 = density_at(dense_r[i - 1]);
        const amrex::Real d1 = density_at(dense_r[i]);
        metric[i] =
            metric[i - 1] + 0.5_rt * (d0 + d1) * (dense_r[i] - dense_r[i - 1]);
    }

    const int n_intervals =
        amrex::max(1, static_cast<int>(std::ceil(metric.back())));
    meta.radius.resize(n_intervals);
    meta.dr.resize(n_intervals);
    meta.chord.resize(n_intervals);
    meta.twist.resize(n_intervals);
    meta.epsilon_profile.resize(n_intervals);
    meta.theta_counts.assign(n_intervals, 1);

    // Place one actuator sample at the center of each equal-metric interval.
    // This gives smaller dr where chord or epsilon requires more radial
    // resolution, while preserving a compact fixed-size velocity-sampling
    // array.
    for (int i = 0; i < n_intervals; ++i) {
        const amrex::Real m0 = metric.back() * static_cast<amrex::Real>(i) /
                               static_cast<amrex::Real>(n_intervals);
        const amrex::Real m1 = metric.back() * static_cast<amrex::Real>(i + 1) /
                               static_cast<amrex::Real>(n_intervals);
        const amrex::Real r0 = interp_metric_to_radius(m0, metric, dense_r);
        const amrex::Real r1 = interp_metric_to_radius(m1, metric, dense_r);
        meta.radius[i] = 0.5_rt * (r0 + r1);
        meta.dr[i] = r1 - r0;
        const amrex::Real mu = (meta.radius[i] - meta.root_radius) /
                               (meta.rotor_radius - meta.root_radius);
        meta.chord[i] =
            ::kynema_sgf::interp::linear(meta.span_locs, meta.chord_inp, mu);
        meta.twist[i] =
            ::kynema_sgf::interp::linear(meta.span_locs, meta.twist_inp, mu);
        meta.epsilon_profile[i] = local_epsilon(meta, meta.chord[i]);
    }
}

vs::Tensor rotation_matrix_from_euler_degrees(const vs::Vector& angles)
{
    return vs::zrot(angles.z()) & vs::yrot(angles.y()) & vs::xrot(angles.x());
}

vs::Tensor orientation_matrix_from_normal(const vs::Vector& normal)
{
    const amrex::Real normal_mag = vs::mag(normal);
    if (normal_mag <= std::numeric_limits<amrex::Real>::epsilon()) {
        amrex::Abort("ActuatorSector: rotor_normal must be nonzero");
    }

    const auto e_normal = normal / normal_mag;
    const auto ref = (std::abs(e_normal.z()) < 0.9_rt) ? vs::Vector::khat()
                                                       : vs::Vector::jhat();
    const auto e_r0 = (ref ^ e_normal).unit();
    const auto e_y = (e_normal ^ e_r0).unit();
    return {e_r0, e_y, e_normal, true};
}

vs::Tensor rotation_matrix_from_vector(const vs::Vector& rotation_vector)
{
    const amrex::Real angle = vs::mag(rotation_vector);
    if (angle <= std::numeric_limits<amrex::Real>::epsilon()) {
        return vs::Tensor::identity();
    }
    return vs::quaternion(
        rotation_vector / angle, -::kynema_sgf::utils::degrees(angle));
}

vs::Tensor orientation_at_time(
    const amrex::Real time,
    const vs::Tensor& initial_orientation,
    const vs::Vector& rotor_angular_velocity)
{
    return rotation_matrix_from_vector(rotor_angular_velocity * time) &
           initial_orientation;
}

void blade_basis(
    const amrex::Real theta,
    const vs::Tensor& orientation,
    vs::Vector& e_r,
    vs::Vector& e_theta,
    vs::Vector& e_normal)
{
    const vs::Vector local_r{std::cos(theta), -std::sin(theta), 0.0_rt};
    const vs::Vector local_t{-std::sin(theta), -std::cos(theta), 0.0_rt};
    const vs::Vector local_n{0.0_rt, 0.0_rt, 1.0_rt};
    e_r = orientation & local_r;
    e_theta = orientation & local_t;
    e_normal = orientation & local_n;
}

vs::Vector blade_position(
    const amrex::Real radius,
    const amrex::Real theta,
    const amrex::Real time,
    const vs::Vector& center0,
    const vs::Vector& translation_velocity,
    const vs::Tensor& initial_orientation,
    const vs::Vector& rotor_angular_velocity)
{
    vs::Vector e_r;
    vs::Vector e_theta;
    vs::Vector e_normal;
    blade_basis(
        theta,
        orientation_at_time(time, initial_orientation, rotor_angular_velocity),
        e_r, e_theta, e_normal);
    return center0 + translation_velocity * time + e_r * radius;
}

void update_midpoint_sample_points(ActuatorSector::DataType& data)
{
    auto& meta = data.meta();
    auto& grid = data.grid();
    const auto& time = data.sim().time();
    // Velocity is sampled at the temporal midpoint, matching the actuator-line
    // convention used elsewhere in the code. The resulting sampled velocity is
    // then used by ComputeForceOp for this step.
    const amrex::Real tmid = 0.5_rt * (time.current_time() + time.new_time());
    const amrex::Real dt = timestep_width(data.sim());
    meta.center = meta.center0 + meta.translation_velocity * tmid;
    meta.azimuth = meta.omega * tmid;
    meta.delta_azimuth = meta.omega * dt;

    const auto orientation = orientation_at_time(
        tmid, orientation_matrix_from_normal(meta.rotor_normal),
        meta.rotor_angular_velocity);

    const int nr = static_cast<int>(meta.radius.size());
    const int nvel = meta.num_blades * nr;
    if (static_cast<int>(grid.vel_pos.size()) != nvel) {
        grid.vel_pos.resize(nvel);
        grid.vel.resize(nvel);
        grid.density.resize(nvel);
    }

    for (int ib = 0; ib < meta.num_blades; ++ib) {
        const amrex::Real phase = ::kynema_sgf::utils::two_pi() *
                                  static_cast<amrex::Real>(ib) /
                                  static_cast<amrex::Real>(meta.num_blades);
        for (int ir = 0; ir < nr; ++ir) {
            const int ip = ib * nr + ir;
            vs::Vector e_r;
            vs::Vector e_theta;
            vs::Vector e_normal;
            blade_basis(
                meta.azimuth + phase, orientation, e_r, e_theta, e_normal);
            grid.vel_pos[ip] = meta.center + e_r * meta.radius[ir];
        }
    }
}

void build_gaussian_table(ActuatorSectorData& meta)
{
    if (meta.gaussian_table_error <= 0.0_rt) {
        amrex::Abort("ActuatorSector: gaussian_table_error must be positive");
    }
    const amrex::Real xmax =
        meta.support_radius_over_epsilon * meta.support_radius_over_epsilon;
    // The table stores exp(-x) for x = (r / epsilon)^2 over the compact support
    // 0 <= x <= support_radius_over_epsilon^2.
    // Linear interpolation error for exp(-x) is bounded by h^2 / 8 because
    // max |d^2 exp(-x) / dx^2| = 1 on x >= 0.
    const amrex::Real max_table_spacing =
        std::sqrt(8.0_rt * meta.gaussian_table_error);
    meta.gaussian_table_nintervals =
        amrex::max(1, static_cast<int>(std::ceil(xmax / max_table_spacing)));
    meta.gaussian_table.resize(meta.gaussian_table_nintervals + 1);
    for (int i = 0; i <= meta.gaussian_table_nintervals; ++i) {
        const amrex::Real x =
            xmax * static_cast<amrex::Real>(i) /
            static_cast<amrex::Real>(meta.gaussian_table_nintervals);
        meta.gaussian_table[i] = std::exp(-x);
    }
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const ActuatorSectorData& meta,
    const ActInfo& info)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) {
        return;
    }

    auto ncf = ncutils::NCFile::create(ncfile, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string nb_name = "num_blades";
    const std::string nr_name = "num_radial_stations";
    const std::vector<std::string> tbr{nt_name, nb_name, nr_name};
    const std::vector<std::string> tbr_vec{nt_name, nb_name, nr_name, "ndim"};

    ncf.enter_def_mode();
    ncf.put_attr("title", "Kynema-SGF actuator-sector blade-load output");
    ncf.put_attr("version", ioutils::kynema_sgf_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);

    auto grp = ncf.def_group(info.label);
    grp.put_attr("force_units", "kinematic force per span: m^3 s^-2");
    grp.put_attr("section_force_units", "kinematic force: m^4 s^-2");
    grp.put_attr(
        "rotor_diameter", std::vector<amrex::Real>{meta.rotor_diameter});
    grp.put_attr(
        "root_radius_fraction",
        std::vector<amrex::Real>{meta.root_radius_fraction});
    grp.put_attr("omega", std::vector<amrex::Real>{meta.omega});
    grp.def_dim(nb_name, static_cast<size_t>(meta.num_blades));
    grp.def_dim(nr_name, meta.radius.size());
    grp.def_var("time", NC_DOUBLE, {nt_name});
    grp.def_var("time_index", NC_INT, {nt_name});
    grp.def_var("thrust", NC_DOUBLE, {nt_name});
    grp.def_var("torque", NC_DOUBLE, {nt_name});
    grp.def_var("radius", NC_DOUBLE, {nr_name});
    grp.def_var("dr", NC_DOUBLE, {nr_name});
    grp.def_var("chord", NC_DOUBLE, {nr_name});
    grp.def_var("twist", NC_DOUBLE, {nr_name});
    grp.def_var("epsilon", NC_DOUBLE, {nr_name});
    grp.def_var("theta_count", NC_INT, {nt_name, nr_name});
    grp.def_var("xyz", NC_DOUBLE, tbr_vec);
    grp.def_var("blade_force", NC_DOUBLE, tbr_vec);
    grp.def_var("section_force", NC_DOUBLE, tbr_vec);
    grp.def_var("vel_rel_theta", NC_DOUBLE, tbr);
    grp.def_var("vel_rel_normal", NC_DOUBLE, tbr);
    grp.def_var("force_theta", NC_DOUBLE, tbr);
    grp.def_var("force_normal", NC_DOUBLE, tbr);
    grp.def_var("aoa", NC_DOUBLE, tbr);
    grp.def_var("cl", NC_DOUBLE, tbr);
    grp.def_var("cd", NC_DOUBLE, tbr);
    grp.var("vel_rel_theta")
        .put_attr("description", "relative velocity along local blade tangent");
    grp.var("vel_rel_normal")
        .put_attr("description", "relative velocity along rotor normal");
    grp.var("vel_rel_theta").put_attr("units", "m s^-1");
    grp.var("vel_rel_normal").put_attr("units", "m s^-1");
    ncf.exit_def_mode();

    const size_t nr = static_cast<size_t>(meta.radius.size());
    grp.var("radius").put(meta.radius.data(), {0}, {nr});
    grp.var("dr").put(meta.dr.data(), {0}, {nr});
    grp.var("chord").put(meta.chord.data(), {0}, {nr});
    grp.var("twist").put(meta.twist.data(), {0}, {nr});
    grp.var("epsilon").put(meta.epsilon_profile.data(), {0}, {nr});
    ncf.close();
#else
    amrex::ignore_unused(ncfile, meta, info);
#endif
}

void write_netcdf(
    const std::string& ncfile,
    const ActuatorSectorData& meta,
    const ActInfo& info,
    const ActGrid& grid,
    const amrex::Real time,
    const int time_index)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) {
        return;
    }

    auto ncf = ncutils::NCFile::open(ncfile, NC_WRITE);
    const size_t nt = ncf.dim("num_time_steps").len();
    const size_t nb = static_cast<size_t>(meta.num_blades);
    const size_t nr = meta.radius.size();
    auto grp = ncf.group(info.label);

    grp.var("time").put(&time, {nt}, {1});
    grp.var("time_index").put(&time_index, {nt}, {1});
    grp.var("thrust").put(&meta.thrust, {nt}, {1});
    grp.var("torque").put(&meta.torque, {nt}, {1});
    grp.var("theta_count").put(meta.theta_counts.data(), {nt, 0}, {1, nr});
    grp.var("xyz").put(
        grid.vel_pos[0].data(), {nt, 0, 0, 0}, {1, nb, nr, AMREX_SPACEDIM});
    grp.var("blade_force")
        .put(
            meta.blade_force[0].data(), {nt, 0, 0, 0},
            {1, nb, nr, AMREX_SPACEDIM});

    VecList section_force(meta.blade_force.size());
    RealList vel_rel_theta(meta.vel_rel.size());
    RealList vel_rel_normal(meta.vel_rel.size());
    for (int ib = 0; ib < meta.num_blades; ++ib) {
        for (int ir = 0; ir < static_cast<int>(nr); ++ir) {
            const int ip = ib * static_cast<int>(nr) + ir;
            section_force[ip] = meta.blade_force[ip] * meta.dr[ir];
            vel_rel_theta[ip] = meta.vel_rel[ip].x();
            vel_rel_normal[ip] = meta.vel_rel[ip].y();
        }
    }
    grp.var("section_force")
        .put(
            section_force[0].data(), {nt, 0, 0, 0},
            {1, nb, nr, AMREX_SPACEDIM});
    grp.var("vel_rel_theta").put(vel_rel_theta.data(), {nt, 0, 0}, {1, nb, nr});
    grp.var("vel_rel_normal")
        .put(vel_rel_normal.data(), {nt, 0, 0}, {1, nb, nr});
    grp.var("force_theta")
        .put(meta.blade_force_theta.data(), {nt, 0, 0}, {1, nb, nr});
    grp.var("force_normal")
        .put(meta.blade_force_normal.data(), {nt, 0, 0}, {1, nb, nr});
    grp.var("aoa").put(meta.aoa.data(), {nt, 0, 0}, {1, nb, nr});
    grp.var("cl").put(meta.cl.data(), {nt, 0, 0}, {1, nb, nr});
    grp.var("cd").put(meta.cd.data(), {nt, 0, 0}, {1, nb, nr});
    ncf.close();
#else
    amrex::ignore_unused(ncfile, meta, info, grid, time, time_index);
#endif
}

} // namespace kynema_sgf::actuator::sector

namespace kynema_sgf::actuator::ops {

void ReadInputsOp<ActuatorSector, ActSrcSector>::operator()(
    ActuatorSector::DataType& data, const utils::ActParser& pp)
{
    auto& meta = data.meta();
    pp.get("rotor_diameter", meta.rotor_diameter);
    pp.get("omega", meta.omega);
    pp.get("airfoil_table", meta.airfoil_file);
    pp.query("airfoil_type", meta.airfoil_type);
    pp.query("num_blades", meta.num_blades);
    pp.query("root_radius_fraction", meta.root_radius_fraction);
    pp.query("center", meta.center0);
    pp.query("translation_velocity", meta.translation_velocity);
    pp.query("rotor_normal", meta.rotor_normal);
    if (pp.contains("rotor_orientation") && !pp.contains("rotor_normal")) {
        amrex::Print()
            << "WARNING: ActuatorSector input 'rotor_orientation' is "
               "deprecated. Use 'rotor_normal' instead.\n";
        vs::Vector rotor_orientation_degrees;
        pp.get("rotor_orientation", rotor_orientation_degrees);
        meta.rotor_normal = sector::rotation_matrix_from_euler_degrees(
                                rotor_orientation_degrees) &
                            vs::Vector::khat();
    }
    pp.query(
        "rotor_rotation_degrees_per_revolution",
        meta.rotor_rotation_degrees_per_revolution);
    if (pp.contains("rotor_rotation_per_revolution")) {
        amrex::Print() << "WARNING: ActuatorSector input "
                          "'rotor_rotation_per_revolution' is deprecated. Use "
                          "'rotor_rotation_degrees_per_revolution' instead.\n";
        pp.get(
            "rotor_rotation_per_revolution",
            meta.rotor_rotation_degrees_per_revolution);
    }
    if (pp.contains("rotor_angular_velocity")) {
        pp.get("rotor_angular_velocity", meta.rotor_angular_velocity);
        meta.user_rotor_angular_velocity = true;
    }
    pp.queryarr("span_locs", meta.span_locs);
    pp.queryarr("chord", meta.chord_inp);
    pp.queryarr("twist", meta.twist_inp);
    if (pp.contains("epsilon")) {
        pp.get("epsilon", meta.epsilon);
        meta.user_epsilon = true;
    }
    pp.query("epsilon_chord", meta.epsilon_chord);
    pp.query("epsilon_min", meta.epsilon_min);
    pp.query("epsilon_dr", meta.epsilon_dr);
    pp.query("epsilon_dl", meta.epsilon_dl);
    pp.query("min_chord_dr", meta.min_chord_dr);
    pp.query("support_radius_over_epsilon", meta.support_radius_over_epsilon);
    pp.query("radial_quadrature", meta.radial_quadrature);
    pp.query("swept_quadrature", meta.swept_quadrature);
    pp.query("gaussian_type", meta.gaussian_type);
    pp.query("gaussian_table_error", meta.gaussian_table_error);

    sector::require_midpoint(meta.radial_quadrature, "radial_quadrature");
    sector::require_midpoint(meta.swept_quadrature, "swept_quadrature");
    if ((amrex::toLower(meta.gaussian_type) != "direct") &&
        (amrex::toLower(meta.gaussian_type) != "table")) {
        amrex::Abort("ActuatorSector: gaussian_type must be direct or table");
    }
    if (!pp.contains("epsilon") && !pp.contains("epsilon_chord")) {
        amrex::Abort(
            "ActuatorSector requires specification of one or both of "
            "'epsilon' or 'epsilon_chord'");
    }
    if (pp.contains("epsilon") && pp.contains("epsilon_chord")) {
        amrex::Print()
            << "WARNING: ActuatorSector inputs 'epsilon' and 'epsilon_chord' "
               "are both specified. Using constant 'epsilon'; "
               "'epsilon_chord' and 'epsilon_min' will be ignored.\n";
    }
    AMREX_ALWAYS_ASSERT(meta.num_blades > 0);
    AMREX_ALWAYS_ASSERT(meta.rotor_diameter > 0.0_rt);
    AMREX_ALWAYS_ASSERT(meta.epsilon >= 0.0_rt);
    AMREX_ALWAYS_ASSERT(meta.epsilon_chord >= 0.0_rt);
    AMREX_ALWAYS_ASSERT(meta.epsilon_min >= 0.0_rt);
    AMREX_ALWAYS_ASSERT(
        (meta.user_epsilon && (meta.epsilon > 0.0_rt)) ||
        (!meta.user_epsilon && (meta.epsilon_chord > 0.0_rt)));
    AMREX_ALWAYS_ASSERT(meta.root_radius_fraction >= 0.0_rt);
    AMREX_ALWAYS_ASSERT(meta.root_radius_fraction < 1.0_rt);
    AMREX_ALWAYS_ASSERT(meta.span_locs.size() == meta.chord_inp.size());
    AMREX_ALWAYS_ASSERT(meta.span_locs.size() == meta.twist_inp.size());

    meta.rotor_radius = 0.5_rt * meta.rotor_diameter;
    meta.root_radius = meta.root_radius_fraction * meta.rotor_radius;
    if (!meta.user_rotor_angular_velocity) {
        meta.rotor_angular_velocity =
            (::kynema_sgf::utils::radians(1.0_rt) * std::abs(meta.omega) /
             ::kynema_sgf::utils::two_pi()) *
            meta.rotor_rotation_degrees_per_revolution;
    }

    const amrex::Real max_eps =
        sector::local_epsilon(meta, sector::max_interp_chord(meta));
    const amrex::Real search_radius =
        meta.rotor_radius + meta.support_radius_over_epsilon * max_eps;
    const auto& geom = data.sim().mesh().Geom(0);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto& c = meta.center0;

    const bool moves = vs::mag(meta.translation_velocity) >
                       std::numeric_limits<amrex::Real>::epsilon();
    const bool rotates = vs::mag(meta.rotor_angular_velocity) >
                         std::numeric_limits<amrex::Real>::epsilon();
    if (moves || rotates) {
        data.info().bound_box = amrex::RealBox(
            moves ? plo[0] : c.x() - search_radius,
            moves ? plo[1] : c.y() - search_radius,
            moves ? plo[2] : c.z() - search_radius,
            moves ? phi[0] : c.x() + search_radius,
            moves ? phi[1] : c.y() + search_radius,
            moves ? phi[2] : c.z() + search_radius);
    } else {
        data.info().bound_box = amrex::RealBox(
            c.x() - search_radius, c.y() - search_radius, c.z() - search_radius,
            c.x() + search_radius, c.y() + search_radius,
            c.z() + search_radius);
    }
}

void InitDataOp<ActuatorSector, ActSrcSector>::operator()(
    ActuatorSector::DataType& data)
{
    auto& meta = data.meta();
    auto& grid = data.grid();
    sector::build_radial_grid(meta);
    if (amrex::toLower(meta.gaussian_type) == "table") {
        sector::build_gaussian_table(meta);
    } else {
        meta.gaussian_table_nintervals = 0;
        meta.gaussian_table.clear();
    }
    meta.aflookup =
        AirfoilLoader::load_airfoil(meta.airfoil_file, meta.airfoil_type);

    const int nr = static_cast<int>(meta.radius.size());
    const int nvel = meta.num_blades * nr;
    grid.resize(nvel, nvel);
    meta.vel_rel.assign(nvel, vs::Vector::zero());
    meta.blade_force.assign(nvel, vs::Vector::zero());
    meta.blade_force_theta.assign(nvel, 0.0_rt);
    meta.blade_force_normal.assign(nvel, 0.0_rt);
    meta.aoa.assign(nvel, 0.0_rt);
    meta.cl.assign(nvel, 0.0_rt);
    meta.cd.assign(nvel, 0.0_rt);
    for (int ip = 0; ip < nvel; ++ip) {
        grid.density[ip] = 1.0_rt;
    }
    sector::update_midpoint_sample_points(data);
}

void UpdatePosOp<ActuatorSector, ActSrcSector>::operator()(
    ActuatorSector::DataType& data)
{
    sector::update_midpoint_sample_points(data);
}

void ComputeForceOp<ActuatorSector, ActSrcSector>::operator()(
    ActuatorSector::DataType& data)
{
    auto& meta = data.meta();
    auto& grid = data.grid();
    const auto& time = data.sim().time();
    const amrex::Real t0 = time.current_time();
    const amrex::Real dt = sector::timestep_width(data.sim());
    const amrex::Real tmid = t0 + 0.5_rt * dt;
    const amrex::Real start_theta = meta.omega * t0;
    const amrex::Real dtheta = meta.omega * dt;
    const auto initial_orientation =
        sector::orientation_matrix_from_normal(meta.rotor_normal);
    const auto mid_orientation = sector::orientation_at_time(
        tmid, initial_orientation, meta.rotor_angular_velocity);

    const int nr = static_cast<int>(meta.radius.size());
    const int nvel = meta.num_blades * nr;
    VecList force_theta_normal(nvel, vs::Vector::zero());
    meta.thrust = 0.0_rt;
    meta.torque = 0.0_rt;

    // First compute blade-section loads at the midpoint sample locations. The
    // force here is per unit span [m^3/s^2] in kinematic units and is converted
    // to section force later by multiplying by dr.
    const amrex::Real spin_sign = (meta.omega >= 0.0_rt) ? 1.0_rt : -1.0_rt;
    for (int ib = 0; ib < meta.num_blades; ++ib) {
        const amrex::Real phase = ::kynema_sgf::utils::two_pi() *
                                  static_cast<amrex::Real>(ib) /
                                  static_cast<amrex::Real>(meta.num_blades);
        for (int ir = 0; ir < nr; ++ir) {
            const int ip = ib * nr + ir;
            vs::Vector e_r;
            vs::Vector e_theta;
            vs::Vector e_normal;
            sector::blade_basis(
                start_theta + 0.5_rt * dtheta + phase, mid_orientation, e_r,
                e_theta, e_normal);

            const auto rel_pos = e_r * meta.radius[ir];
            const auto frame_vel = meta.rotor_angular_velocity ^ rel_pos;
            const auto spin_vel = e_theta * (meta.omega * meta.radius[ir]);
            const auto blade_vel =
                meta.translation_velocity + frame_vel + spin_vel;
            const auto rel_wind = grid.vel[ip] - blade_vel;
            const amrex::Real vtheta = rel_wind & e_theta;
            const amrex::Real vnormal = rel_wind & e_normal;
            const auto vplane = e_theta * vtheta + e_normal * vnormal;
            const amrex::Real vmag = vs::mag(vplane);
            // AoA is measured in the local tangential-normal blade plane using
            // a rotation-aware tangential direction. This keeps a mirrored
            // propeller (negative omega and negative twist) aerodynamically
            // equivalent to the positive-rotation propeller, while preserving
            // opposite tangential force for torque cancellation.
            const amrex::Real raw_aoa =
                std::atan2(vnormal, -spin_sign * vtheta) +
                spin_sign * ::kynema_sgf::utils::radians(meta.twist[ir]);
            const amrex::Real aoa =
                std::remainder(raw_aoa, ::kynema_sgf::utils::two_pi());

            amrex::Real cl = 0.0_rt;
            amrex::Real cd = 0.0_rt;
            (*meta.aflookup)(aoa, cl, cd);

            const amrex::Real qval = 0.5_rt * vmag * vmag * meta.chord[ir];
            const amrex::Real lift = qval * cl;
            const amrex::Real drag = qval * cd;
            const auto drag_dir =
                (vmag > std::numeric_limits<amrex::Real>::epsilon())
                    ? vplane.unit()
                    : e_theta;
            const auto lift_dir = spin_sign * (drag_dir ^ e_r).unit();
            const auto force_on_fluid =
                -((lift_dir * lift) + (drag_dir * drag));
            const amrex::Real ftheta = force_on_fluid & e_theta;
            const amrex::Real fnormal = force_on_fluid & e_normal;
            force_theta_normal[ip] = {ftheta, fnormal, 0.0_rt};

            meta.blade_force[ip] = force_on_fluid;
            meta.blade_force_theta[ip] = ftheta;
            meta.blade_force_normal[ip] = fnormal;
            meta.vel_rel[ip] = {vtheta, vnormal, 0.0_rt};
            meta.aoa[ip] = ::kynema_sgf::utils::degrees(aoa);
            meta.cl[ip] = cl;
            meta.cd[ip] = cd;
            meta.thrust += fnormal * meta.dr[ir];
            meta.torque += ftheta * meta.radius[ir] * meta.dr[ir];
        }
    }

    // Convert each radial blade-section load into one or more force projection
    // points over the swept sector for this timestep. If dt is small enough,
    // theta_counts is one and the sector naturally reduces to actuator-line
    // behavior.
    int nforce = 0;
    for (int ir = 0; ir < nr; ++ir) {
        const amrex::Real swept_speed =
            vs::mag(meta.translation_velocity) +
            std::abs(meta.omega) * meta.radius[ir] +
            vs::mag(meta.rotor_angular_velocity) * meta.radius[ir];
        const int ntheta = amrex::max(
            1, static_cast<int>(std::ceil(
                   swept_speed * dt * meta.epsilon_dl /
                   meta.epsilon_profile[ir])));
        meta.theta_counts[ir] = ntheta;
        nforce += meta.num_blades * ntheta;
    }

    grid.pos.resize(nforce);
    grid.force.resize(nforce);
    grid.epsilon.resize(nforce);
    grid.orientation.assign(nforce, vs::Tensor::identity());
    int iq = 0;
    for (int ib = 0; ib < meta.num_blades; ++ib) {
        const amrex::Real phase = ::kynema_sgf::utils::two_pi() *
                                  static_cast<amrex::Real>(ib) /
                                  static_cast<amrex::Real>(meta.num_blades);
        for (int ir = 0; ir < nr; ++ir) {
            const int ip = ib * nr + ir;
            const int ntheta = meta.theta_counts[ir];
            for (int it = 0; it < ntheta; ++it) {
                const amrex::Real xi = (static_cast<amrex::Real>(it) + 0.5_rt) /
                                       static_cast<amrex::Real>(ntheta);
                const amrex::Real t = t0 + xi * dt;
                const amrex::Real theta = start_theta + xi * dtheta + phase;
                const auto orientation = sector::orientation_at_time(
                    t, initial_orientation, meta.rotor_angular_velocity);
                vs::Vector e_r;
                vs::Vector e_theta;
                vs::Vector e_normal;
                sector::blade_basis(theta, orientation, e_r, e_theta, e_normal);
                grid.pos[iq] = meta.center0 + meta.translation_velocity * t +
                               e_r * meta.radius[ir];
                // Split the radial section force uniformly across the swept
                // quadrature points so the integrated force remains unchanged.
                const amrex::Real wt =
                    meta.dr[ir] / static_cast<amrex::Real>(ntheta);
                grid.force[iq] = wt * ((e_theta * force_theta_normal[ip].x()) +
                                       (e_normal * force_theta_normal[ip].y()));
                grid.epsilon[iq] = vs::Vector::one() * meta.epsilon_profile[ir];
                ++iq;
            }
        }
    }
}

void ProcessOutputsOp<ActuatorSector, ActSrcSector>::read_io_options(
    const utils::ActParser& pp)
{
    pp.query("output_frequency", m_out_freq);
}

void ProcessOutputsOp<ActuatorSector, ActSrcSector>::prepare_outputs(
    const std::string& out_dir)
{
    m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
    sector::prepare_netcdf_file(m_nc_filename, m_data.meta(), m_data.info());
}

void ProcessOutputsOp<ActuatorSector, ActSrcSector>::write_outputs()
{
    const auto& time = m_data.sim().time();
    const int tidx = time.time_index();
    if ((m_out_freq <= 0) || (tidx % m_out_freq != 0)) {
        return;
    }

    sector::write_netcdf(
        m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
        time.new_time(), tidx);
}

void ActSrcOp<ActuatorSector, ActSrcSector>::initialize() { copy_to_device(); }

void ActSrcOp<ActuatorSector, ActSrcSector>::copy_to_device()
{
    const auto& grid = m_data.grid();
    const int npts = static_cast<int>(grid.pos.size());
    if (static_cast<int>(m_pos.size()) != npts) {
        m_pos.resize(npts);
        m_force.resize(npts);
        m_epsilon.resize(npts);
    }
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
        m_pos.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.force.begin(), grid.force.end(),
        m_force.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.epsilon.begin(), grid.epsilon.end(),
        m_epsilon.begin());
    const auto& table = m_data.meta().gaussian_table;
    if (static_cast<int>(m_gaussian_table.size()) !=
        static_cast<int>(table.size())) {
        m_gaussian_table.resize(table.size());
    }
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, table.begin(), table.end(),
        m_gaussian_table.begin());
}

void ActSrcOp<ActuatorSector, ActSrcSector>::operator()(
    const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    BL_PROFILE("kynema-sgf::ActSrcOp<ActuatorSector>");

    const auto& bx = mfi.tilebox();
    const auto bxa =
        kynema_sgf::utils::realbox_to_box(m_data.info().bound_box, geom);
    const auto& bxi = bx & bxa;
    if (bxi.isEmpty()) {
        return;
    }

    const auto& sarr = m_act_src(lev).array(mfi);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const int npts = static_cast<int>(m_pos.size());
    const auto* pos = m_pos.data();
    const auto* force = m_force.data();
    const auto* eps = m_epsilon.data();

    const auto& meta = m_data.meta();
    const amrex::Real support2 =
        meta.support_radius_over_epsilon * meta.support_radius_over_epsilon;
    const bool use_table = amrex::toLower(meta.gaussian_type) == "table";
    const auto* table = m_gaussian_table.data();
    const int table_nintervals = meta.gaussian_table_nintervals;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const vs::Vector cc{
            problo[0] + ((i + 0.5_rt) * dx[0]),
            problo[1] + ((j + 0.5_rt) * dx[1]),
            problo[2] + ((k + 0.5_rt) * dx[2]),
        };

        amrex::RealArray src_force = {0.0_rt};
        for (int ip = 0; ip < npts; ++ip) {
            const auto dist = cc - pos[ip];
            const vs::Vector rr{
                dist.x() / eps[ip].x(), dist.y() / eps[ip].y(),
                dist.z() / eps[ip].z()};
            const amrex::Real rr_sqr = vs::mag_sqr(rr);
            if (rr_sqr >= support2) {
                continue;
            }

            amrex::Real eval = 0.0_rt;
            if (use_table) {
                // rr_sqr already equals x = (r / epsilon)^2. The table is
                // uniform in x over [0, support_radius_over_epsilon^2].
                const amrex::Real scaled =
                    rr_sqr / support2 *
                    static_cast<amrex::Real>(table_nintervals);
                const int ilo = static_cast<int>(scaled);
                const int ihi = amrex::min(ilo + 1, table_nintervals);
                const amrex::Real frac = scaled - static_cast<amrex::Real>(ilo);
                eval = (1.0_rt - frac) * table[ilo] + frac * table[ihi];
            } else {
                eval = std::exp(-rr_sqr);
            }
            const amrex::Real eps_fac = eps[ip].x() * eps[ip].y() * eps[ip].z();
            // Normalized 3D Gaussian. The eps product gives the physical
            // volume scaling, and gaussian_norm_3d() makes the infinite-domain
            // kernel integrate to one.
            const amrex::Real gauss_fac =
                sector::gaussian_norm_3d() / eps_fac * eval;
            src_force[0] += gauss_fac * force[ip].x();
            src_force[1] += gauss_fac * force[ip].y();
            src_force[2] += gauss_fac * force[ip].z();
        }

        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}

} // namespace kynema_sgf::actuator::ops
