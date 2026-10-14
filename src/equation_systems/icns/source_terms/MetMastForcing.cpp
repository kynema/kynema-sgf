#include "src/equation_systems/icns/source_terms/MetMastForcing.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

#include <fstream>
#include <limits>

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

namespace {

template <typename T>
void copy_to_device(
    const amrex::Vector<T>& host, amrex::Gpu::DeviceVector<T>& device)
{
    device.resize(host.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, host.begin(), host.end(), device.begin());
}

} // namespace

// Named namespace: nvcc rejects device lambdas in functions with internal
// linkage
namespace metmast {

//! Device view of the station profiles
struct StationView
{
    const amrex::Real* x;
    const amrex::Real* y;
    const int* offset;
    const amrex::Real* z;
    const amrex::Real* vel;
    const amrex::Real* sigma;
    int nstations;
};

struct ForcingParams
{
    amrex::Real inv_rh2;
    amrex::Real inv_rz2;
    amrex::Real cutoff;
    amrex::Real sigma_factor;
    amrex::Real inv_tau;
    amrex::Real max_rate;
};

/** Add the met-mast forcing to the momentum source
 *
 *  With terrain, the height is measured from the local terrain and the cells
 *  inside the terrain are not forced.
 */
template <bool HasTerrain>
void add_forcing(
    const amrex::Geometry& geom,
    const StationView& st,
    const ForcingParams& prm,
    const amrex::MultiArray4<amrex::Real const>& vel_arrs,
    const amrex::MultiArray4<amrex::Real const>& terrain_arrs,
    const amrex::MultiArray4<int const>& blank_arrs,
    amrex::MultiFab& src_term)
{
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            if constexpr (HasTerrain) {
                if (blank_arrs[nbx](i, j, k) == 1) {
                    return;
                }
                z = amrex::max<amrex::Real>(
                    z - terrain_arrs[nbx](i, j, k), 0.5_rt * dx[2]);
            }
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const auto& vel = vel_arrs[nbx];

            amrex::Real sum_w = 0.0_rt;
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> sum_wt{
                AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)};
            for (int s = 0; s < st.nstations; ++s) {
                const amrex::Real xs = x - st.x[s];
                const amrex::Real ys = y - st.y[s];
                const amrex::Real rh2 = ((xs * xs) + (ys * ys)) * prm.inv_rh2;
                if (rh2 > prm.cutoff) {
                    continue;
                }
                // Bracket the height in the station profile; outside the
                // measured range use the end value and taper the weight
                const int lbot = st.offset[s];
                const int ltop = st.offset[s + 1] - 1;
                int la = lbot;
                int lb = lbot;
                amrex::Real frac = 0.0_rt;
                amrex::Real dz_out = 0.0_rt;
                if (z <= st.z[lbot]) {
                    dz_out = st.z[lbot] - z;
                } else if (z >= st.z[ltop]) {
                    la = ltop;
                    lb = ltop;
                    dz_out = z - st.z[ltop];
                } else {
                    lb = lbot + 1;
                    while (st.z[lb] < z) {
                        ++lb;
                    }
                    la = lb - 1;
                    frac = (z - st.z[la]) / (st.z[lb] - st.z[la]);
                }
                const amrex::Real ri2 = rh2 + (dz_out * dz_out * prm.inv_rz2);
                if (ri2 > prm.cutoff) {
                    continue;
                }
                const amrex::Real w = std::exp(-0.25_rt * ri2);
                sum_w += w;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    const int ia = (la * AMREX_SPACEDIM) + n;
                    const int ib = (lb * AMREX_SPACEDIM) + n;
                    const amrex::Real ref =
                        ((1.0_rt - frac) * st.vel[ia]) + (frac * st.vel[ib]);
                    const amrex::Real band =
                        prm.sigma_factor * (((1.0_rt - frac) * st.sigma[ia]) +
                                            (frac * st.sigma[ib]));
                    const amrex::Real target = amrex::min(
                        amrex::max(vel(i, j, k, n), ref - band), ref + band);
                    sum_wt[n] += w * target;
                }
            }
            if (sum_w <= 0.0_rt) {
                return;
            }
            const amrex::Real rate = amrex::min(
                amrex::min(sum_w, 1.0_rt) * prm.inv_tau, prm.max_rate);
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                src_arrs[nbx](i, j, k, n) -=
                    rate * (vel(i, j, k, n) - (sum_wt[n] / sum_w));
            }
        });
}

} // namespace metmast

MetMastForcing::MetMastForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_sim(sim)
{
    amrex::ParmParse pp_abl("ABL");
    std::string point_file;
    amrex::Vector<std::string> profile_files;
    pp_abl.query("metmast_1dprofile_file", point_file);
    pp_abl.queryarr("metmast_profile_files", profile_files);
    if (point_file.empty() && profile_files.empty()) {
        amrex::Abort(
            "MetMastForcing: set ABL.metmast_1dprofile_file and/or "
            "ABL.metmast_profile_files");
    }
    if (!point_file.empty()) {
        read_point_file(point_file);
    }
    for (const auto& fname : profile_files) {
        read_profile_file(fname);
    }
    if (num_stations() == 0) {
        amrex::Abort(
            "MetMastForcing: no measurements found in the input files");
    }

    pp_abl.query("meso_timescale", m_timescale);
    pp_abl.query("metmast_timescale", m_timescale);
    pp_abl.query("metmast_horizontal_radius", m_horizontal_radius);
    pp_abl.query("metmast_vertical_radius", m_vertical_radius);
    pp_abl.query("metmast_damping_radius", m_damping_radius);
    pp_abl.query("metmast_sigma_factor", m_sigma_factor);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_timescale > 0.0_rt, "MetMastForcing: time scale must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (m_horizontal_radius > 0.0_rt) && (m_vertical_radius > 0.0_rt),
        "MetMastForcing: radii must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_sigma_factor >= 0.0_rt,
        "MetMastForcing: sigma factor must not be negative");

    copy_to_device(m_station_x, m_station_x_d);
    copy_to_device(m_station_y, m_station_y_d);
    copy_to_device(m_level_offset, m_level_offset_d);
    copy_to_device(m_level_z, m_level_z_d);
    copy_to_device(m_level_vel, m_level_vel_d);
    copy_to_device(m_level_sigma, m_level_sigma_d);
}

MetMastForcing::~MetMastForcing() = default;

void MetMastForcing::add_level(
    const amrex::Real z,
    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& vel,
    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& sigma)
{
    m_level_z.push_back(z);
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        m_level_vel.push_back(vel[n]);
        m_level_sigma.push_back(sigma[n]);
    }
}

void MetMastForcing::read_point_file(const std::string& fname)
{
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("MetMastForcing: cannot find met-mast file " + fname);
    }
    //! x y z u v w T; the temperature is not used
    amrex::Real x, y, z, u, v, w, temp;
    while (ifh >> x >> y >> z >> u >> v >> w >> temp) {
        m_station_x.push_back(x);
        m_station_y.push_back(y);
        add_level(z, {u, v, w}, {0.0_rt, 0.0_rt, 0.0_rt});
        m_level_offset.push_back(static_cast<int>(m_level_z.size()));
    }
    if (!ifh.eof()) {
        amrex::Abort("MetMastForcing: cannot parse met-mast file " + fname);
    }
}

void MetMastForcing::read_profile_file(const std::string& fname)
{
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("MetMastForcing: cannot find profile file " + fname);
    }
    amrex::Real x, y;
    if (!(ifh >> x >> y)) {
        amrex::Abort("MetMastForcing: cannot read location in " + fname);
    }
    //! z u v w su sv sw, with z above the terrain in ascending order
    const auto nlevels_old = m_level_z.size();
    amrex::Real z, u, v, w, su, sv, sw;
    while (ifh >> z >> u >> v >> w >> su >> sv >> sw) {
        if ((m_level_z.size() > nlevels_old) && (z <= m_level_z.back())) {
            amrex::Abort(
                "MetMastForcing: heights must be increasing in " + fname);
        }
        if ((su < 0.0_rt) || (sv < 0.0_rt) || (sw < 0.0_rt)) {
            amrex::Abort(
                "MetMastForcing: negative standard deviation in " + fname);
        }
        add_level(z, {u, v, w}, {su, sv, sw});
    }
    if (!ifh.eof()) {
        amrex::Abort("MetMastForcing: cannot parse profile file " + fname);
    }
    if (m_level_z.size() == nlevels_old) {
        amrex::Abort("MetMastForcing: no heights in profile file " + fname);
    }
    m_station_x.push_back(x);
    m_station_y.push_back(y);
    m_level_offset.push_back(static_cast<int>(m_level_z.size()));
}

void MetMastForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const metmast::StationView st{
        .x = m_station_x_d.data(),
        .y = m_station_y_d.data(),
        .offset = m_level_offset_d.data(),
        .z = m_level_z_d.data(),
        .vel = m_level_vel_d.data(),
        .sigma = m_level_sigma_d.data(),
        .nstations = num_stations()};

    const amrex::Real dt = m_time.delta_t();
    const metmast::ForcingParams prm{
        .inv_rh2 = 1.0_rt / (m_horizontal_radius * m_horizontal_radius),
        .inv_rz2 = 1.0_rt / (m_vertical_radius * m_vertical_radius),
        .cutoff = m_damping_radius,
        .sigma_factor = m_sigma_factor,
        .inv_tau = 1.0_rt / m_timescale,
        .max_rate = (dt > 0.0_rt) ? 1.0_rt / dt
                                  : std::numeric_limits<amrex::Real>::max()};

    const auto& geom = m_mesh.Geom(lev);
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();

    const auto& repo = m_sim.repo();
    if (repo.field_exists("terrain_height")) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            repo.int_field_exists("terrain_blank"),
            "MetMastForcing: terrain_height requires terrain_blank");
        metmast::add_forcing<true>(
            geom, st, prm, vel_arrs,
            repo.get_field("terrain_height")(lev).const_arrays(),
            repo.get_int_field("terrain_blank")(lev).const_arrays(), src_term);
    } else {
        metmast::add_forcing<false>(
            geom, st, prm, vel_arrs, amrex::MultiArray4<amrex::Real const>(),
            amrex::MultiArray4<int const>(), src_term);
    }
}

} // namespace kynema_sgf::pde::icns
