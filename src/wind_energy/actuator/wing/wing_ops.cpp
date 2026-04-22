#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::wing {

void init_data_structures(WingBaseData& wdata, ActGrid& grid)
{
    int npts = wdata.num_pts;
    grid.resize(npts);

    // Wing span
    auto wspan = wdata.end - wdata.start;
    // Compute transformation matrix as a quaternion rotation about wing span
    auto tmat = vs::quaternion(wspan, wdata.pitch);
    // Epsilon (chord, thickness, span) -> (chord, span, thickness)
    const auto& epsin = wdata.eps_inp;
    vs::Vector eps{epsin.x(), epsin.z(), epsin.y()};

    // Equal spacing along span
    auto dx = (1.0_rt / static_cast<amrex::Real>(npts - 1)) * wspan;

    for (int i = 0; i < npts; ++i) {
        grid.pos[i] = wdata.start + static_cast<amrex::Real>(i) * dx;
    }

    // Initialize remaining data
    grid.epsilon.assign(npts, eps);
    grid.orientation.assign(npts, tmat);
    grid.force.assign(npts, vs::Vector::zero());
    grid.vel_pos.assign(grid.pos.begin(), grid.pos.end());
    grid.vel.assign(npts, vs::Vector::zero());

    // Assign length of actuator segments
    wdata.dx.assign(npts, vs::mag(dx));
    // The first and last segments have half width
    wdata.dx.front() *= 0.5_rt;
    wdata.dx.back() *= 0.5_rt;

    wdata.vel_rel.assign(npts, vs::Vector::zero());
    wdata.aoa.assign(npts, 0.0_rt);
    wdata.cl.assign(npts, 0.0_rt);
    wdata.cd.assign(npts, 0.0_rt);
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const WingBaseData& meta,
    const ActInfo& info,
    const ActGrid& grid)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) {
        return;
    }

    auto ncf = ncutils::NCFile::create(ncfile, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_actuator_points";
    const std::vector<std::string> two_dim{nt_name, np_name};

    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind fixed wing actuator output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);

    auto grp = ncf.def_group(info.label);
    // clang-format off
    grp.put_attr("epsilon",
        std::vector<amrex::Real>{meta.eps_inp.x(),
                         meta.eps_inp.y(), meta.eps_inp.z()});
    // clang-format on
    grp.def_dim(np_name, meta.num_pts);
    grp.def_var("time", NC_DOUBLE, {nt_name});
    grp.def_var("pitch", NC_DOUBLE, {nt_name});
    grp.def_var("integrated_lift", NC_DOUBLE, {nt_name});
    grp.def_var("integrated_drag", NC_DOUBLE, {nt_name});
    grp.def_var("xyz", NC_DOUBLE, {np_name, "ndim"});
    grp.def_var("chord", NC_DOUBLE, {np_name});
    grp.def_var("dx", NC_DOUBLE, {np_name});
    grp.def_var("veff", NC_DOUBLE, {nt_name, np_name, "ndim"});
    grp.def_var("vrel", NC_DOUBLE, {nt_name, np_name, "ndim"});
    grp.def_var("body_force", NC_DOUBLE, {nt_name, np_name, "ndim"});
    grp.def_var("aoa", NC_DOUBLE, two_dim);
    grp.def_var("cl", NC_DOUBLE, two_dim);
    grp.def_var("cd", NC_DOUBLE, two_dim);
    ncf.exit_def_mode();

    {
        const auto npts = static_cast<size_t>(meta.num_pts);
        const std::vector<size_t> start{0, 0};
        const std::vector<size_t> count{npts, AMREX_SPACEDIM};
        auto xyz = grp.var("xyz");
        xyz.put(grid.pos[0].data(), start, count);
        auto chord = grp.var("chord");
        chord.put(meta.chord.data(), {0}, {npts});
        auto dx = grp.var("dx");
        dx.put(meta.dx.data(), {0}, {npts});
    }
#else
    amrex::ignore_unused(ncfile, meta, info, grid);
#endif
}

void write_netcdf(
    const std::string& ncfile,
    const WingBaseData& meta,
    const ActInfo& info,
    const ActGrid& grid,
    const amrex::Real time)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) {
        return;
    }

    auto ncf = ncutils::NCFile::open(ncfile, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of next timestep
    const size_t nt = ncf.dim(nt_name).len();
    const auto npts = static_cast<size_t>(meta.num_pts);

    std::vector<size_t> start{nt, 0};
    std::vector<size_t> count{1, npts};
    auto grp = ncf.group(info.label);
    grp.var("time").put(&time, {nt}, {1});
    grp.var("pitch").put(&meta.pitch, {nt}, {1});
    grp.var("integrated_lift").put(&meta.lift, {nt}, {1});
    grp.var("integrated_drag").put(&meta.drag, {nt}, {1});
    grp.var("vrel").put(
        meta.vel_rel[0].data(), {nt, 0, 0}, {1, npts, AMREX_SPACEDIM});
    grp.var("veff").put(
        grid.vel[0].data(), {nt, 0, 0}, {1, npts, AMREX_SPACEDIM});
    grp.var("body_force")
        .put(grid.force[0].data(), {nt, 0, 0}, {1, npts, AMREX_SPACEDIM});
    grp.var("aoa").put(meta.aoa.data(), start, count);
    grp.var("cl").put(meta.cl.data(), start, count);
    grp.var("cd").put(meta.cd.data(), start, count);
#else
    amrex::ignore_unused(ncfile, meta, info, grid, time);
#endif
}

void refresh_wing_position(VecList& vpoints, VecList fpoints, const int npts)
{
    for (int ip = 0; ip < npts; ++ip) {
        // Move velocity points to latest force points
        vpoints[ip].x() = fpoints[ip].x();
        vpoints[ip].y() = fpoints[ip].y();
        vpoints[ip].z() = fpoints[ip].z();
    }
}

void new_wing_position_velocity(
    VecList& points,
    vs::Vector& vtr,
    const int npts,
    const amrex::Real tn,
    const amrex::Real tnp1,
    const std::string& motion,
    const amrex::Real period,
    const vs::Vector svec)
{
    // Get displacement of points from n to n+1
    // Also, if translation velocity changes, update it
    vs::Vector disp{0.0_rt, 0.0_rt, 0.0_rt};
    // Do nothing for "none"
    if (amrex::toLower(motion) == "linear") {
        // Use velocity to get displacement
        disp.x() = vtr.x() * (tnp1 - tn);
        disp.y() = vtr.y() * (tnp1 - tn);
        disp.z() = vtr.z() * (tnp1 - tn);
        // Velocity is unchanged
    } else if (amrex::toLower(motion) == "sine") {
        // Calculate displacement using sine
        disp.x() =
            svec.x() *
            (std::sin(
                 2.0_rt * std::numbers::pi_v<amrex::Real> * tnp1 / period) -
             std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * tn / period));
        disp.y() =
            svec.y() *
            (std::sin(
                 2.0_rt * std::numbers::pi_v<amrex::Real> * tnp1 / period) -
             std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * tn / period));
        disp.z() =
            svec.z() *
            (std::sin(
                 2.0_rt * std::numbers::pi_v<amrex::Real> * tnp1 / period) -
             std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * tn / period));
        // The translational velocity over the time step is disp/dt
        vtr.x() = disp.x() /
                  (tnp1 - tn + std::numeric_limits<amrex::Real>::epsilon());
        vtr.y() = disp.y() /
                  (tnp1 - tn + std::numeric_limits<amrex::Real>::epsilon());
        vtr.z() = disp.z() /
                  (tnp1 - tn + std::numeric_limits<amrex::Real>::epsilon());
        // The tiny number in the denominator is important for initialization
    }
    for (int ip = 0; ip < npts; ++ip) {
        // Move points according to displacement
        points[ip].x() = points[ip].x() + disp.x();
        points[ip].y() = points[ip].y() + disp.y();
        points[ip].z() = points[ip].z() + disp.z();
    }
}

} // namespace amr_wind::actuator::wing
