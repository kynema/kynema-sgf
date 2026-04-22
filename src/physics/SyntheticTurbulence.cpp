#include <memory>

#include "amr-wind/physics/SyntheticTurbulence.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
namespace synth_turb {

struct LinearShearOp
{
    const amrex::Real m_hmin;
    const amrex::Real m_hmax;
    const amrex::Real m_vstart;
    const amrex::Real m_vstop;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real ht) const
    {
        amrex::Real vel = m_vstart + ((m_vstop - m_vstart) * (ht - m_hmin) /
                                      (m_hmax - m_hmin));
        return amrex::max(m_vstart, amrex::min(vel, m_vstop));
    }
};

class LinearShearProfile : public MeanProfile
{
private:
public:
    LinearShearProfile(
        amrex::Real h_min,
        amrex::Real h_max,
        amrex::Real vel_start,
        amrex::Real vel_stop,
        int shear_dir)
        : MeanProfile(0.5_rt * (vel_start + vel_stop), shear_dir)
        , m_op{
              .m_hmin = h_min,
              .m_hmax = h_max,
              .m_vstart = vel_start,
              .m_vstop = vel_stop}
    {}

    ~LinearShearProfile() override = default;

    // cppcheck-suppress duplInheritedMember
    [[nodiscard]] LinearShearOp device_instance() const { return m_op; }

private:
    LinearShearOp m_op;
};

struct PowerLawOp
{
    const amrex::Real m_ref_vel;
    const amrex::Real m_ref_height;
    const amrex::Real m_alpha;
    const amrex::Real m_hoffset;
    const amrex::Real m_umin;
    const amrex::Real m_umax;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    operator()(amrex::Real height) const
    {
        const amrex::Real heff = height - m_hoffset;
        amrex::Real pfac =
            (heff > 0.0_rt) ? std::pow((heff / m_ref_height), m_alpha) : 0.0_rt;
        return m_ref_vel * amrex::min(amrex::max(pfac, m_umin), m_umax);
    }
};

class PowerLawProfile : public MeanProfile
{
public:
    PowerLawProfile(
        amrex::Real ref_vel,
        amrex::Real ref_height,
        amrex::Real alpha,
        int shear_dir,
        amrex::Real h_offset,
        amrex::Real umin,
        amrex::Real umax)
        : MeanProfile(ref_vel, shear_dir)
        , m_op{
              .m_ref_vel = ref_vel,
              .m_ref_height = ref_height,
              .m_alpha = alpha,
              .m_hoffset = h_offset,
              .m_umin = umin,
              .m_umax = umax}
    {}

    ~PowerLawProfile() override = default;

    // cppcheck-suppress duplInheritedMember
    [[nodiscard]] PowerLawOp device_instance() const { return m_op; }

private:
    PowerLawOp m_op;
};

} // namespace synth_turb

namespace {

/** Parse the NetCDF turbulence database and determine details of the turbulence
 *box.
 *
 *. Initializes the dimensions and grid length, sizes in SynthTurbData. Also
 *  allocates the necessary memory for the perturbation velocities.
 *
 *. \param turbFile Information regarding NetCDF data identifiers
 *. \param turbGrid Turbulence data
 */
void process_nc_file(const std::string& turb_filename, SynthTurbData& turb_grid)
{
#ifdef AMR_WIND_USE_NETCDF
    auto ncf = ncutils::NCFile::open(turb_filename, NC_NOWRITE);

    // Grid dimensions
    AMREX_ASSERT(ncf.dim("ndim").len() == AMREX_SPACEDIM);
    auto nx = ncf.dim("nx").len();
    auto ny = ncf.dim("ny").len();
    auto nz = ncf.dim("nz").len();

    turb_grid.box_dims[0] = static_cast<int>(nx);
    turb_grid.box_dims[1] = static_cast<int>(ny);
    turb_grid.box_dims[2] = static_cast<int>(nz);

    // Box lengths and resolution
    auto box_len = ncf.var("box_lengths");
    box_len.get(turb_grid.box_len.data());
    auto dx = ncf.var("dx");
    dx.get(turb_grid.dx.data());

    ncf.close();

    // Create data structures to store the perturbation velocities for two
    // planes
    const size_t grid_size = 2 * ny * nz;
    turb_grid.uvel.resize(grid_size);
    turb_grid.vvel.resize(grid_size);
    turb_grid.wvel.resize(grid_size);
    turb_grid.uvel_d.resize(grid_size);
    turb_grid.vvel_d.resize(grid_size);
    turb_grid.wvel_d.resize(grid_size);
#else
    amrex::ignore_unused(turb_filename, turb_grid);
#endif
}

/** Load two planes of data that bound the current timestep
 *
 *  The data for the y and z directions are loaded for the entire grid at the
 * two planes
 */
void load_turb_plane_data(
    const std::string& turb_filename,
    SynthTurbData& turb_grid,
    const int il,
    const int ir)
{
    BL_PROFILE("amr-wind::SyntheticTurbulence::load_plane_data");
#ifdef AMR_WIND_USE_NETCDF
    auto ncf = ncutils::NCFile::open(turb_filename, NC_NOWRITE);

    // clang-format off
    std::vector<size_t> start{static_cast<size_t>(il), 0, 0};
    std::vector<size_t> count{2, static_cast<size_t>(turb_grid.box_dims[1]),
                               static_cast<size_t>(turb_grid.box_dims[2])};
    // clang-format on

    auto uvel = ncf.var("uvel");
    auto vvel = ncf.var("vvel");
    auto wvel = ncf.var("wvel");

    if ((ir - il) == 1) {
        // two consecutive planes load them in one shot
        uvel.get(turb_grid.uvel.data(), start, count);
        vvel.get(turb_grid.vvel.data(), start, count);
        wvel.get(turb_grid.wvel.data(), start, count);
    } else {
        // Load the planes separately
        count[0] = 1;
        uvel.get(turb_grid.uvel.data(), start, count);
        vvel.get(turb_grid.vvel.data(), start, count);
        wvel.get(turb_grid.wvel.data(), start, count);

        start[0] = static_cast<size_t>(ir);
        const size_t offset = static_cast<size_t>(turb_grid.box_dims[1]) *
                              static_cast<size_t>(turb_grid.box_dims[2]);
        uvel.get(&turb_grid.uvel[offset], start, count);
        vvel.get(&turb_grid.vvel[offset], start, count);
        wvel.get(&turb_grid.wvel[offset], start, count);
    }

    // Update left and right indices for future checks
    turb_grid.ileft = il;
    turb_grid.iright = ir;

    ncf.close();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, turb_grid.uvel.begin(), turb_grid.uvel.end(),
        turb_grid.uvel_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, turb_grid.vvel.begin(), turb_grid.vvel.end(),
        turb_grid.vvel_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, turb_grid.wvel.begin(), turb_grid.wvel.end(),
        turb_grid.wvel_d.begin());
#else
    amrex::ignore_unused(turb_filename, turb_grid, il, ir);
#endif
}

/** Determine the left/right indices for a given point along a particular
 * direction
 *
 *  \param turb_grid Turbulence box data
 *  \param dir Direction of search (0 = x, 1 = y, 2 = z)
 *  \param xin Coordinate value in local coordinate frame corresponding to
 * direction provided
 *  \param il Index of the lower bound (populated by this function)
 *  \param ir Index of the upper bound (populated by this function)
 */
void get_lr_indices(
    const SynthTurbData& turb_grid,
    const int dir,
    const amrex::Real xin,
    int& il,
    int& ir)
{
    const amrex::Real xbox = xin - (std::floor(xin / turb_grid.box_len[dir]) *
                                    turb_grid.box_len[dir]);

    il = static_cast<int>(std::floor(xbox / turb_grid.dx[dir]));
    ir = il + 1;
    if (ir >= turb_grid.box_dims[dir]) {
        ir -= turb_grid.box_dims[dir];
    }
}

/** Determine the left/right indices for a given point along a particular
 * direction
 *
 *  This overload also populates the fractions of left/right states to be used
 *  for interpolations.
 *
 *  \param turb_grid Turbulence box data
 *  \param dir Direction of search (0 = x, 1 = y, 2 = z)
 *  \param xin Coordinate value in local coordinate frame corresponding to
 *  direction provided
 *  \param il Index of the lower bound (populated by this function)
 *  \param ir Index of the upper bound (populated by this function)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_lr_indices(
    const SynthTurbDeviceData& turb_grid,
    const int dir,
    const amrex::Real xin,
    int& il,
    int& ir,
    amrex::Real& rxl,
    amrex::Real& rxr)
{
    const amrex::Real xbox = xin - (std::floor(xin / turb_grid.box_len[dir]) *
                                    turb_grid.box_len[dir]);

    il = static_cast<int>(std::floor(xbox / turb_grid.dx[dir]));
    ir = il + 1;
    if (ir >= turb_grid.box_dims[dir]) {
        ir -= turb_grid.box_dims[dir];
    }

    const amrex::Real xfrac = xbox - (turb_grid.dx[dir] * il);
    rxr = xfrac / turb_grid.dx[dir];
    rxl = (1.0_rt - rxr);
}

/** Determine if a given point (in local frame) is within the turbulence box
 *
 *  If the point is found within the box, also determine the indices and
 *  interpolation weights for the y and z directions.
 *
 *  @return True if the point is inside the 2-D box
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool find_point_in_box(
    const SynthTurbDeviceData& t_grid, const vs::Vector& pt, InterpWeights& wt)
{
    // Get y and z w.r.t. the lower corner of the grid
    const amrex::Real yy = pt[1] + (t_grid.box_len[1] * 0.5_rt);
    const amrex::Real zz = pt[2] + (t_grid.box_len[2] * 0.5_rt);
    bool inBox =
        ((yy >= 0.0_rt) && (yy <= t_grid.box_len[1]) && (zz >= 0.0_rt) &&
         (zz <= t_grid.box_len[2]));

    if (inBox) {
        get_lr_indices(t_grid, 1, yy, wt.jl, wt.jr, wt.yl, wt.yr);
        get_lr_indices(t_grid, 2, zz, wt.kl, wt.kr, wt.zl, wt.zr);
    }

    return inBox;
}

/** Interpolate the perturbation velocity to a given point from the grid data
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void interp_perturb_vel(
    const SynthTurbDeviceData& t_grid, const InterpWeights& wt, vs::Vector& vel)
{
    const int nz = t_grid.box_dims[2];
    const int nynz = t_grid.box_dims[1] * t_grid.box_dims[2];
    // Indices of the 2-D cell that contains the sampling point
    amrex::Array<int, 4> qidx = {
        (wt.jl * nz) + wt.kl, (wt.jr * nz) + wt.kl, (wt.jr * nz) + wt.kr,
        (wt.jl * nz) + wt.kr};

    vs::Vector vel_l, vel_r;

    // Left quad (t = t)
    vel_l[0] = (wt.yl * wt.zl * t_grid.uvel[qidx[0]]) +
               (wt.yr * wt.zl * t_grid.uvel[qidx[1]]) +
               (wt.yr * wt.zr * t_grid.uvel[qidx[2]]) +
               (wt.yl * wt.zr * t_grid.uvel[qidx[3]]);
    vel_l[1] = (wt.yl * wt.zl * t_grid.vvel[qidx[0]]) +
               (wt.yr * wt.zl * t_grid.vvel[qidx[1]]) +
               (wt.yr * wt.zr * t_grid.vvel[qidx[2]]) +
               (wt.yl * wt.zr * t_grid.vvel[qidx[3]]);
    vel_l[2] = (wt.yl * wt.zl * t_grid.wvel[qidx[0]]) +
               (wt.yr * wt.zl * t_grid.wvel[qidx[1]]) +
               (wt.yr * wt.zr * t_grid.wvel[qidx[2]]) +
               (wt.yl * wt.zr * t_grid.wvel[qidx[3]]);

    for (int& i : qidx) {
        i += nynz;
    }

    // Right quad (t = t+delta_t)
    vel_r[0] = (wt.yl * wt.zl * t_grid.uvel[qidx[0]]) +
               (wt.yr * wt.zl * t_grid.uvel[qidx[1]]) +
               (wt.yr * wt.zr * t_grid.uvel[qidx[2]]) +
               (wt.yl * wt.zr * t_grid.uvel[qidx[3]]);
    vel_r[1] = (wt.yl * wt.zl * t_grid.vvel[qidx[0]]) +
               (wt.yr * wt.zl * t_grid.vvel[qidx[1]]) +
               (wt.yr * wt.zr * t_grid.vvel[qidx[2]]) +
               (wt.yl * wt.zr * t_grid.vvel[qidx[3]]);
    vel_r[2] = (wt.yl * wt.zl * t_grid.wvel[qidx[0]]) +
               (wt.yr * wt.zl * t_grid.wvel[qidx[1]]) +
               (wt.yr * wt.zr * t_grid.wvel[qidx[2]]) +
               (wt.yl * wt.zr * t_grid.wvel[qidx[3]]);

    // Interpolation in time
    vel = wt.xl * vel_l + wt.xr * vel_r;
}

} // namespace

SyntheticTurbulence::SyntheticTurbulence(const CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_turb_force(sim.repo().declare_field("synth_turb_forcing", 3))
{
#ifndef AMR_WIND_USE_NETCDF
    amrex::Abort(
        "SyntheticTurbulence: AMR-Wind was not built with NetCDF support.");
#endif
    amrex::ParmParse pp("SynthTurb");

    // NetCDF file containing the turbulence data
    pp.query("turbulence_file", m_turb_filename);
    process_nc_file(m_turb_filename, m_turb_grid);

    // Load position and orientation of the grid
    amrex::Real wind_direction{270.};
    pp.query("wind_direction", wind_direction);
    amrex::Vector<amrex::Real> location{0.0_rt, 0.0_rt, 0.0_rt};
    pp.queryarr("grid_location", location);

    std::string mean_wind_type = "ConstValue";
    pp.get("mean_wind_type", mean_wind_type);

    if (mean_wind_type == "ConstValue") {
        amrex::ParmParse pp_vel("ConstValue.velocity");
        amrex::Vector<amrex::Real> vel;
        pp_vel.getarr("value", vel);
        amrex::Real wind_speed = vs::mag(vs::Vector{vel[0], vel[1], vel[2]});
        m_wind_profile =
            std::make_unique<synth_turb::MeanProfile>(wind_speed, 2);
    } else if (mean_wind_type == "LinearProfile") {
        amrex::ParmParse pp_vel("LinearProfile.velocity");

        amrex::Real zmin, zmax;
        pp_vel.get("start", zmin);
        pp_vel.get("stop", zmax);

        amrex::Vector<amrex::Real> start_val, stop_val;
        pp_vel.getarr("start_val", start_val);
        pp_vel.getarr("stop_val", stop_val);

        int shear_dir = 2;
        pp_vel.query("direction", shear_dir);

        m_wind_profile = std::make_unique<synth_turb::LinearShearProfile>(
            zmin, zmax,
            vs::mag(vs::Vector{start_val[0], start_val[1], start_val[2]}),
            vs::mag(vs::Vector{stop_val[0], stop_val[1], stop_val[2]}),
            shear_dir);
    } else if (mean_wind_type == "PowerLawProfile") {
        amrex::ParmParse pp_vel("PowerLawProfile.velocity");

        // Default reference height is the center of the turbulence grid
        amrex::Real zref = location[2];
        pp.get("zref", zref);
        amrex::Vector<amrex::Real> vel;
        pp.getarr("uref", vel);
        amrex::Real wind_speed = vs::mag(vs::Vector{vel[0], vel[1], vel[2]});

        amrex::Real alpha, zoffset, umin, umax;
        pp.get("shear_exponent", alpha);
        pp.get("zoffset", zoffset);
        pp.get("umin", umin);
        pp.get("umax", umax);

        int shear_dir = 2;
        pp_vel.query("direction", shear_dir);

        umin /= wind_speed;
        umax /= wind_speed;
        m_wind_profile = std::make_unique<synth_turb::PowerLawProfile>(
            wind_speed, zref, alpha, shear_dir, zoffset, umin, umax);
    } else {
        amrex::Abort(
            "SyntheticTurbulence: invalid mean wind type specified = " +
            mean_wind_type);
    }

    m_mean_wind_type = mean_wind_type;
    // Smearing factors
    pp.get("gauss_smearing_factor", m_epsilon);
    m_gauss_scaling =
        1.0_rt / (m_epsilon * std::sqrt(std::numbers::pi_v<amrex::Real>));

    // Time offsets if any...
    pp.query("time_offset", m_time_offset);

    // Duration
    pp.query("duration", m_duration);

    // Done reading user inputs, process derived data

    // Center of the grid
    m_turb_grid.origin[0] = location[0];
    m_turb_grid.origin[1] = location[1];
    m_turb_grid.origin[2] = location[2];

    // Compute box-fixed reference frame.
    //
    // x-direction points to flow direction (convert from compass direction to
    // vector)
    m_turb_grid.tr_mat = vs::zrot(270.0_rt - wind_direction);

    amrex::Print() << "Synthetic turbulence forcing initialized \n"
                   << "  Turbulence file = " << m_turb_filename << "\n"
                   << "  Box lengths = [" << m_turb_grid.box_len[0] << ", "
                   << m_turb_grid.box_len[1] << ", " << m_turb_grid.box_len[2]
                   << "]\n"
                   << "  Box dims = [" << m_turb_grid.box_dims[0] << ", "
                   << m_turb_grid.box_dims[1] << ", " << m_turb_grid.box_dims[2]
                   << "]\n"
                   << "  Grid dx = [" << m_turb_grid.dx[0] << ", "
                   << m_turb_grid.dx[1] << ", " << m_turb_grid.dx[2] << "]\n"
                   << "  Centroid (forcing plane) = [" << m_turb_grid.origin[0]
                   << ", " << m_turb_grid.origin[1] << ", "
                   << m_turb_grid.origin[2] << "]\n"
                   << "  Mean wind profile: U = "
                   << m_wind_profile->reference_velocity()
                   << " m/s; Dir = " << wind_direction
                   << " deg; type = " << mean_wind_type << '\n';
}

void SyntheticTurbulence::initialize_fields(
    int /*level*/, const amrex::Geometry& /*geom*/)
{}

void SyntheticTurbulence::pre_advance_work()
{
    if (m_is_init) {
        initialize();
    }

    update();
}

void SyntheticTurbulence::initialize()
{
    BL_PROFILE("amr-wind::SyntheticTurbulence::initialize");
    // Convert current time to an equivalent length based on the reference
    // velocity to determine the position within the turbulence grid
    const amrex::Real curTime = m_time.new_time() - m_time_offset;
    const amrex::Real eqivLen = m_wind_profile->reference_velocity() * curTime;
    int il, ir;
    get_lr_indices(m_turb_grid, 0, eqivLen, il, ir);
    load_turb_plane_data(m_turb_filename, m_turb_grid, il, ir);

    m_is_init = false;
}

void SyntheticTurbulence::update()
{
    BL_PROFILE("amr-wind::SyntheticTurbulence::update");

    // Convert current time to an equivalent length based on the reference
    // velocity to determine the position within the turbulence grid
    const amrex::Real cur_time = m_time.new_time() - m_time_offset;
    const amrex::Real eqiv_len =
        m_wind_profile->reference_velocity() * cur_time;

    // Stop update if the current time is past the requested duration of the
    // injection of the synthetic turbulence and if this duration is positive
    if (m_duration > 0.0_rt && cur_time > m_duration) {
        const amrex::Vector<amrex::Real> zeros{0.0_rt, 0.0_rt, 0.0_rt};
        m_turb_force.setVal(zeros, m_turb_force.num_grow()[0]);
        return;
    }

    InterpWeights weights;
    SynthTurbDeviceData turb_grid(m_turb_grid);
    get_lr_indices(
        turb_grid, 0, eqiv_len, weights.il, weights.ir, weights.xl, weights.xr);

    // Check if we need to refresh the planes
    if (weights.il != m_turb_grid.ileft) {
        load_turb_plane_data(
            m_turb_filename, m_turb_grid, weights.il, weights.ir);
    }

    if (m_mean_wind_type == "ConstValue") {
        update_impl(turb_grid, weights, m_wind_profile->device_instance());
    } else if (m_mean_wind_type == "LinearProfile") {
        const auto* vfunc =
            dynamic_cast<synth_turb::LinearShearProfile*>(m_wind_profile.get());
        if (vfunc != nullptr) {
            update_impl(turb_grid, weights, vfunc->device_instance());
        }
    } else if (m_mean_wind_type == "PowerLawProfile") {
        const auto* vfunc =
            dynamic_cast<synth_turb::PowerLawProfile*>(m_wind_profile.get());
        if (vfunc != nullptr) {
            update_impl(turb_grid, weights, vfunc->device_instance());
        }
    }
}

template <typename T>
void SyntheticTurbulence::update_impl(
    const SynthTurbDeviceData& turb_grid,
    const InterpWeights& weights,
    const T& velfunc)
{
    const auto& repo = m_turb_force.repo();
    const auto& geom_vec = repo.mesh().Geom();

    const int sdir = (*m_wind_profile).shear_dir();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = geom.ProbLoArray();

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const auto& trmat = m_turb_grid.tr_mat;
        const auto& origin = m_turb_grid.origin;
        const auto& gauss_scaling = m_gauss_scaling;
        const auto& epsilon = m_epsilon;

        const auto& turb_force_arrs = m_turb_force(lev).arrays();
        const auto& rho_arrs = m_density(lev).arrays();

        amrex::ParallelFor(
            m_turb_force(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // Position vector in local turbulence grid frame
                vs::Vector xyz_l;
                // velocity in local frame
                vs::Vector vel_l;
                // velocity in global frame
                vs::Vector vel_g;

                vs::Vector xyz_g{
                    problo[0] + ((i + 0.5_rt) * dx),
                    problo[1] + ((j + 0.5_rt) * dy),
                    problo[2] + ((k + 0.5_rt) * dz)};

                // Transform position vector from global inertial
                // reference frame to local reference frame attached to
                // the turbulence grid.
                xyz_l = trmat & (xyz_g - origin);

                InterpWeights wts_loc = weights;

                // Check if the point is in the box, if not we skip this
                // node. The function will also populate the interpolation
                // weights for points that are determined to be within the
                // box.
                const bool ptInBox =
                    find_point_in_box(turb_grid, xyz_l, wts_loc);
                if (ptInBox) {
                    // Interpolate perturbation velocities in the local
                    // reference frame
                    interp_perturb_vel(turb_grid, wts_loc, vel_l);
                    // Transform velocity vector from local reference
                    // frame back to the global inertial frame
                    vel_g = vel_l & trmat;

                    // Based on the equations in
                    // http://doi.wiley.com/10.1002/we.1608
                    // v_n in Eq. 10
                    const amrex::Real v_mag = vs::mag(vel_g);
                    // (V_n + 1/2 v_n) in Eq. 10
                    const amrex::Real v_mag_total =
                        (velfunc(xyz_g[sdir]) + (0.5_rt * v_mag));

                    // Smearing factor (see Eq. 11). The normal direction to
                    // the grid is the x-axis of the local reference frame
                    // by construction
                    const amrex::Real term1 = xyz_l[0] / epsilon;
                    const amrex::Real eta =
                        std::exp(-(term1 * term1)) * gauss_scaling;
                    const amrex::Real factor = v_mag_total * eta;

                    turb_force_arrs[nbx](i, j, k, 0) =
                        rho_arrs[nbx](i, j, k) * vel_g[0] * factor;
                    turb_force_arrs[nbx](i, j, k, 1) =
                        rho_arrs[nbx](i, j, k) * vel_g[1] * factor;
                    turb_force_arrs[nbx](i, j, k, 2) =
                        rho_arrs[nbx](i, j, k) * vel_g[2] * factor;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
