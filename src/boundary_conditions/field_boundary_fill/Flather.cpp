#include "src/CFDSim.H"
#include "src/boundary_conditions/field_boundary_fill/Flather.H"
#include "src/boundary_conditions/field_boundary_fill/FillFlather.H"
#include "src/utilities/index_operations.H"
#include "src/utilities/constants.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_GpuAtomic.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"
#include <algorithm>
#include <cmath>

using namespace amrex::literals;

namespace kynema_sgf {

Flather::Flather(CFDSim& sim)
    : m_sim(sim)
    , m_time(m_sim.time())
    , m_repo(m_sim.repo())
    , m_mesh(m_sim.mesh())
    , m_velocity(m_sim.repo().get_field("velocity"))
    , m_u_mac(m_sim.repo().get_field("u_mac"))
    , m_v_mac(m_sim.repo().get_field("v_mac"))
    , m_vof(m_sim.repo().get_field("vof"))
{
    if (!m_repo.field_exists("vof")) {
        amrex::Abort("Flather BC requires the vof field");
    }
    if (m_repo.int_field_exists("terrain_blank")) {
        m_terrain_blank = &m_repo.get_int_field("terrain_blank");
    }

    {
        amrex::ParmParse pp(identifier());
        pp.query("max_velocity_scale_factor", m_velocity_scale_factor_max);
        pp.query("min_velocity_scale_factor", m_velocity_scale_factor_min);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    // Vector at the x boundary extends in y direction
    // Vector at the y boundary extends in x direction
    m_xlo_uhliq.resize(1, m_mesh.Geom());
    m_xhi_uhliq.resize(1, m_mesh.Geom());
    m_ylo_uhliq.resize(0, m_mesh.Geom());
    m_yhi_uhliq.resize(0, m_mesh.Geom());
    m_xlo_uhmix.resize(1, m_mesh.Geom());
    m_xhi_uhmix.resize(1, m_mesh.Geom());
    m_ylo_uhmix.resize(0, m_mesh.Geom());
    m_yhi_uhmix.resize(0, m_mesh.Geom());
    m_xlo_bnd_uh.resize(1, m_mesh.Geom());
    m_xhi_bnd_uh.resize(1, m_mesh.Geom());
    m_ylo_bnd_uh.resize(0, m_mesh.Geom());
    m_yhi_bnd_uh.resize(0, m_mesh.Geom());

    m_xlo_h.resize(1, m_mesh.Geom());
    m_xhi_h.resize(1, m_mesh.Geom());
    m_ylo_h.resize(0, m_mesh.Geom());
    m_yhi_h.resize(0, m_mesh.Geom());
    m_xlo_bnd_h.resize(1, m_mesh.Geom());
    m_xhi_bnd_h.resize(1, m_mesh.Geom());
    m_ylo_bnd_h.resize(0, m_mesh.Geom());
    m_yhi_bnd_h.resize(0, m_mesh.Geom());
}

void Flather::post_init_actions()
{
    m_velocity.add_fill_patch_op<FillFlather>(m_mesh, m_time, *this);
    compute_internal_z_averages();
}

void Flather::pre_advance_work() { compute_internal_z_averages(); }

void Flather::accumulate_boundary(
    int current_level,
    int idir,
    int phase_switch,
    bool is_low,
    MultiLevelVector& out_uvec,
    MultiLevelVector& out_hvec,
    bool sample_boundary,
    FieldState fstate,
    bool use_mac_fields) const
{
    AMREX_ALWAYS_ASSERT(idir == 0 || idir == 1);

    auto& int_h = out_uvec.host_data(current_level);
    auto& dist_h = out_hvec.host_data(current_level);
    const int nline = static_cast<int>(int_h.size());

    amrex::Gpu::DeviceVector<amrex::Real> uh_sum_d(nline, 0.0_rt);
    amrex::Gpu::DeviceVector<amrex::Real> vof_sum_d(nline, 0.0_rt);

    const amrex::Real tiny = constants::TIGHT_TOL;

    // Loop through current level and all below
    for (int lev = current_level; lev >= 0; --lev) {

        // Cumulative refinement ratio, tangent to the boundary, between this
        // level and the current level
        int rr = 1;
        for (int flev = lev; flev < current_level; ++flev) {
            rr *= m_mesh.refRatio(flev)[1 - idir];
        }

        amrex::iMultiFab level_mask;
        if (lev < current_level) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& geom = m_mesh.Geom(lev);
        const auto dz = geom.CellSizeArray()[2];
        const auto& dom = geom.Domain();
        const int bidx = is_low ? dom.smallEnd(idir) : dom.bigEnd(idir);
        const int shift_to_boundary = sample_boundary ? (is_low ? -1 : 1) : 0;
        const auto& src_vel =
            use_mac_fields ? ((idir == 0) ? m_u_mac : m_v_mac) : m_velocity;
        if (shift_to_boundary != 0) {
            AMREX_ALWAYS_ASSERT(src_vel.num_grow()[idir] > 0);
            AMREX_ALWAYS_ASSERT(m_vof.num_grow()[idir] > 0);
        }

        const auto& vel_mf =
            src_vel.state(use_mac_fields ? FieldState::New : fstate)(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_mf, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            // Iterate in cell-centered index space, even for MAC fields, so
            // the box type matches the domain and the ii_v/jj_v face offsets
            amrex::Box bx = mfi.tilebox(amrex::IntVect::TheZeroVector()) & dom;
            if (!bx.ok()) {
                continue;
            }
            if (bidx < bx.smallEnd(idir) || bidx > bx.bigEnd(idir)) {
                continue;
            }

            // Limit box to along boundary
            bx.setSmall(idir, bidx);
            bx.setBig(idir, bidx);

            const auto vel_arr = vel_mf.const_array(mfi);
            // Due to the order that fillphysbc is called in prepare_boundaries
            // (velocity, then scalars), the vof data is not available at
            // alternate (nph) states when this is called.
            const auto vof_arr = m_vof(lev).const_array(mfi);
            const auto mask_arr = level_mask.const_array(mfi);
            const bool use_terrain = (m_terrain_blank != nullptr);
            const auto terrain_blank_arr =
                use_terrain ? (*m_terrain_blank)(lev).const_array(mfi)
                            : amrex::Array4<int const>();

            amrex::Real* uh_sum = uh_sum_d.data();
            amrex::Real* vof_sum = vof_sum_d.data();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const int ii = (idir == 0) ? (i + shift_to_boundary) : i;
                const int jj = (idir == 1) ? (j + shift_to_boundary) : j;
                const int ii_v =
                    ii +
                    static_cast<int>(use_mac_fields && idir == 0 && !is_low);
                const int jj_v =
                    jj +
                    static_cast<int>(use_mac_fields && idir == 1 && !is_low);

                if (use_terrain && terrain_blank_arr(i, j, k) != 0) {
                    return;
                }
                // This index is tangent to the boundary
                const int idx_lev = (idir == 0) ? j : i;
                // Convert to current level indices
                const int idx_min = idx_lev * rr;
                const int idx_max = idx_min + rr - 1;

                for (int idx = idx_min; idx <= idx_max; ++idx) {
                    const auto liquid_height =
                        vof_arr(ii, jj, k) * dz * mask_arr(i, j, k);
                    // Threads share these sums, so the atomic must be
                    // host-safe as well as device-safe
                    amrex::HostDevice::Atomic::Add(
                        &vof_sum[idx], liquid_height);
                    auto vel_height = liquid_height;
                    if (phase_switch == 0 &&
                        vof_arr(ii, jj, k) < 1.0_rt - tiny) {
                        // Needs to be fully liquid to be counted
                        vel_height = 0.0_rt;
                    }
                    if (phase_switch == 1 &&
                        vof_arr(ii, jj, k) >= 1.0_rt - tiny) {
                        // Needs to be mixture to be counted
                        vel_height = 0.0_rt;
                    }
                    // Fully gas cells never count because of vof multiplier
                    amrex::Real local_vel =
                        vel_arr(ii_v, jj_v, k, use_mac_fields ? 0 : idir);
                    // If summing internal velocities, only allow outflow
                    if (!sample_boundary) {
                        local_vel = is_low ? amrex::min(local_vel, 0.0_rt)
                                           : amrex::max(local_vel, 0.0_rt);
                    }
                    amrex::HostDevice::Atomic::Add(
                        &uh_sum[idx], local_vel * vel_height);
                }
            });
        }
    }

    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, uh_sum_d.begin(), uh_sum_d.end(),
        int_h.begin());
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, vof_sum_d.begin(), vof_sum_d.end(),
        dist_h.begin());

    amrex::ParallelDescriptor::ReduceRealSum(int_h.data(), nline);
    amrex::ParallelDescriptor::ReduceRealSum(dist_h.data(), nline);
}

void Flather::compute_internal_z_averages()
{
    BL_PROFILE("kynema-sgf::Flather::compute_internal_z_averages");

    const int nlevels = m_repo.num_active_levels();
    const int nlevels_geom = static_cast<int>(m_mesh.Geom().size());
    if (m_xlo_uhliq.size() != nlevels_geom) {
        m_xlo_uhliq.resize(1, m_mesh.Geom());
        m_xhi_uhliq.resize(1, m_mesh.Geom());
        m_ylo_uhliq.resize(0, m_mesh.Geom());
        m_yhi_uhliq.resize(0, m_mesh.Geom());
        m_xlo_uhmix.resize(1, m_mesh.Geom());
        m_xhi_uhmix.resize(1, m_mesh.Geom());
        m_ylo_uhmix.resize(0, m_mesh.Geom());
        m_yhi_uhmix.resize(0, m_mesh.Geom());
        m_xlo_bnd_uh.resize(1, m_mesh.Geom());
        m_xhi_bnd_uh.resize(1, m_mesh.Geom());
        m_ylo_bnd_uh.resize(0, m_mesh.Geom());
        m_yhi_bnd_uh.resize(0, m_mesh.Geom());

        m_xlo_h.resize(1, m_mesh.Geom());
        m_xhi_h.resize(1, m_mesh.Geom());
        m_ylo_h.resize(0, m_mesh.Geom());
        m_yhi_h.resize(0, m_mesh.Geom());
        m_xlo_bnd_h.resize(1, m_mesh.Geom());
        m_xhi_bnd_h.resize(1, m_mesh.Geom());
        m_ylo_bnd_h.resize(0, m_mesh.Geom());
        m_yhi_bnd_h.resize(0, m_mesh.Geom());
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        this->accumulate_boundary(
            lev, 0, 0, true, m_xlo_uhliq, m_xlo_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 0, 0, false, m_xhi_uhliq, m_xhi_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 1, 0, true, m_ylo_uhliq, m_ylo_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 1, 0, false, m_yhi_uhliq, m_yhi_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 0, 1, true, m_xlo_uhmix, m_xlo_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 0, 1, false, m_xhi_uhmix, m_xhi_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 1, 1, true, m_ylo_uhmix, m_ylo_h, false, FieldState::New);
        this->accumulate_boundary(
            lev, 1, 1, false, m_yhi_uhmix, m_yhi_h, false, FieldState::New);
    }

    m_xlo_uhliq.copy_host_to_device();
    m_xhi_uhliq.copy_host_to_device();
    m_ylo_uhliq.copy_host_to_device();
    m_yhi_uhliq.copy_host_to_device();

    m_xlo_uhmix.copy_host_to_device();
    m_xhi_uhmix.copy_host_to_device();
    m_ylo_uhmix.copy_host_to_device();
    m_yhi_uhmix.copy_host_to_device();

    m_xlo_h.copy_host_to_device();
    m_xhi_h.copy_host_to_device();
    m_ylo_h.copy_host_to_device();
    m_yhi_h.copy_host_to_device();
    amrex::Gpu::streamSynchronize();
}

void Flather::compute_boundary_z_averages(
    int lev, FieldState fstate, bool use_mac_fields)
{
    BL_PROFILE("kynema-sgf::Flather::compute_boundary_z_averages");

    // accumulating boundaries needs to happen prior to applying fillpatch op

    this->accumulate_boundary(
        lev, 0, -1, true, m_xlo_bnd_uh, m_xlo_bnd_h, true, fstate,
        use_mac_fields);
    this->accumulate_boundary(
        lev, 0, -1, false, m_xhi_bnd_uh, m_xhi_bnd_h, true, fstate,
        use_mac_fields);
    this->accumulate_boundary(
        lev, 1, -1, true, m_ylo_bnd_uh, m_ylo_bnd_h, true, fstate,
        use_mac_fields);
    this->accumulate_boundary(
        lev, 1, -1, false, m_yhi_bnd_uh, m_yhi_bnd_h, true, fstate,
        use_mac_fields);

    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_xlo_bnd_uh.host_data(lev).begin(),
        m_xlo_bnd_uh.host_data(lev).end(),
        m_xlo_bnd_uh.device_data(lev).begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_xhi_bnd_uh.host_data(lev).begin(),
        m_xhi_bnd_uh.host_data(lev).end(),
        m_xhi_bnd_uh.device_data(lev).begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_ylo_bnd_uh.host_data(lev).begin(),
        m_ylo_bnd_uh.host_data(lev).end(),
        m_ylo_bnd_uh.device_data(lev).begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_yhi_bnd_uh.host_data(lev).begin(),
        m_yhi_bnd_uh.host_data(lev).end(),
        m_yhi_bnd_uh.device_data(lev).begin());

    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_xlo_bnd_h.host_data(lev).begin(),
        m_xlo_bnd_h.host_data(lev).end(), m_xlo_bnd_h.device_data(lev).begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_xhi_bnd_h.host_data(lev).begin(),
        m_xhi_bnd_h.host_data(lev).end(), m_xhi_bnd_h.device_data(lev).begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_ylo_bnd_h.host_data(lev).begin(),
        m_ylo_bnd_h.host_data(lev).end(), m_ylo_bnd_h.device_data(lev).begin());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_yhi_bnd_h.host_data(lev).begin(),
        m_yhi_bnd_h.host_data(lev).end(), m_yhi_bnd_h.device_data(lev).begin());
    amrex::Gpu::streamSynchronize();
}

void Flather::set_velocity(
    const int lev,
    const amrex::Real time,
    const Field& fld,
    amrex::MultiFab& mfab,
    const int /* dcomp */,
    const int orig_comp) const
{
    BL_PROFILE("kynema-sgf::Flather::set_velocity");

    const auto& geom = m_mesh.Geom(lev);
    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const int numcomp = mfab.nComp();
    const auto& domain = geom.growPeriodicDomain(nghost);

    const auto fstate = time > 0.0_rt ? FieldState::Old : FieldState::New;
    const amrex::Real tiny = constants::TIGHT_TOL;
    const amrex::Real v_threshold = 1.0e-6_rt;
    const auto grav_z = -m_gravity[2];
    const auto vscale_min = m_velocity_scale_factor_min;
    const auto vscale_max = m_velocity_scale_factor_max;

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow)) {
            continue;
        }

        const int idir = ori.coordDir();
        if (idir == 2) {
            continue;
        }
        // Check if orientation aligns with supplied field and choose component
        bool skip_fill = false;
        int fcomp = 0;
        if (numcomp == 1) {
            // MAC velocity: only normal field is valid
            skip_fill = (orig_comp != idir);
            // Only a single component available
        } else {
            // Cell-centered velocity: field always valid
            fcomp = idir;
            // Select valid component
        }
        if (skip_fill) {
            continue;
        }

        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);
        const auto shift_to_interior =
            amrex::IntVect::TheDimensionVector(idir) * (ori.isLow() ? 1 : -1);
        const amrex::Real* xlo_uhliq = m_xlo_uhliq.device_data(lev).data();
        const amrex::Real* xhi_uhliq = m_xhi_uhliq.device_data(lev).data();
        const amrex::Real* ylo_uhliq = m_ylo_uhliq.device_data(lev).data();
        const amrex::Real* yhi_uhliq = m_yhi_uhliq.device_data(lev).data();
        const amrex::Real* xlo_uhmix = m_xlo_uhmix.device_data(lev).data();
        const amrex::Real* xhi_uhmix = m_xhi_uhmix.device_data(lev).data();
        const amrex::Real* ylo_uhmix = m_ylo_uhmix.device_data(lev).data();
        const amrex::Real* yhi_uhmix = m_yhi_uhmix.device_data(lev).data();
        const amrex::Real* xlo_h = m_xlo_h.device_data(lev).data();
        const amrex::Real* xhi_h = m_xhi_h.device_data(lev).data();
        const amrex::Real* ylo_h = m_ylo_h.device_data(lev).data();
        const amrex::Real* yhi_h = m_yhi_h.device_data(lev).data();
        const amrex::Real* xlo_bnd_uh = m_xlo_bnd_uh.device_data(lev).data();
        const amrex::Real* xhi_bnd_uh = m_xhi_bnd_uh.device_data(lev).data();
        const amrex::Real* ylo_bnd_uh = m_ylo_bnd_uh.device_data(lev).data();
        const amrex::Real* yhi_bnd_uh = m_yhi_bnd_uh.device_data(lev).data();
        const amrex::Real* xlo_bnd_h = m_xlo_bnd_h.device_data(lev).data();
        const amrex::Real* xhi_bnd_h = m_xhi_bnd_h.device_data(lev).data();
        const amrex::Real* ylo_bnd_h = m_ylo_bnd_h.device_data(lev).data();
        const amrex::Real* yhi_bnd_h = m_yhi_bnd_h.device_data(lev).data();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            auto shift_to_cc = amrex::IntVect(0);
            const auto& bx = utils::face_aware_boundary_box_intersection(
                shift_to_cc, gbx, dbx, ori);
            if (!bx.ok()) {
                continue;
            }

            const auto& arr = mfab[mfi].array();
            const auto& ref_arr = fld.state(fstate)(lev)[mfi].const_array();

            const auto vof_arr = m_vof.state(fstate)(lev).const_array(mfi);
            const bool use_terrain = (m_terrain_blank != nullptr);
            const auto terrain_blank_arr =
                use_terrain ? (*m_terrain_blank)(lev).const_array(mfi)
                            : amrex::Array4<int const>();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::IntVect iv{i, j, k};
                const amrex::IntVect iv_cc = iv + shift_to_cc;
                const amrex::IntVect iv_adj = iv + shift_to_interior;
                const amrex::IntVect iv_adj_cc = iv_adj + shift_to_cc;

                amrex::Real boundary_val = arr(iv, fcomp);
                amrex::Real interior_liq = arr(iv_adj, fcomp);
                amrex::Real interior_mix = arr(iv_adj, fcomp);
                amrex::Real boundary_h = 0.0_rt;
                amrex::Real interior_h = 0.0_rt;
                // Vectors at x boundaries extend in y direction
                // Vectors at y boundaries extend in x direction
                if (idir == 0) {
                    interior_liq = ori.isLow() ? xlo_uhliq[iv_adj[1]]
                                               : xhi_uhliq[iv_adj[1]];
                    interior_mix = ori.isLow() ? xlo_uhmix[iv_adj[1]]
                                               : xhi_uhmix[iv_adj[1]];
                    interior_h =
                        ori.isLow() ? xlo_h[iv_adj[1]] : xhi_h[iv_adj[1]];
                    boundary_val =
                        ori.isLow() ? xlo_bnd_uh[iv[1]] : xhi_bnd_uh[iv[1]];
                    boundary_h =
                        ori.isLow() ? xlo_bnd_h[iv[1]] : xhi_bnd_h[iv[1]];
                } else {
                    interior_liq = ori.isLow() ? ylo_uhliq[iv_adj[0]]
                                               : yhi_uhliq[iv_adj[0]];
                    interior_mix = ori.isLow() ? ylo_uhmix[iv_adj[0]]
                                               : yhi_uhmix[iv_adj[0]];
                    interior_h =
                        ori.isLow() ? ylo_h[iv_adj[0]] : yhi_h[iv_adj[0]];
                    boundary_val =
                        ori.isLow() ? ylo_bnd_uh[iv[0]] : yhi_bnd_uh[iv[0]];
                    boundary_h =
                        ori.isLow() ? ylo_bnd_h[iv[0]] : yhi_bnd_h[iv[0]];
                }

                // Set velocity to zero and skip for terrain
                if (use_terrain && terrain_blank_arr(iv_adj_cc) != 0) {
                    arr(iv, fcomp) = 0.0_rt;
                    return;
                }
                // Check interior or boundary vof for liquid
                const amrex::Real interior_vof = vof_arr(iv_adj_cc);
                const amrex::Real boundary_vof = vof_arr(iv_cc);
                // Do nothing here if not liquid
                if ((interior_vof < tiny) && (boundary_vof < tiny)) {
                    return;
                }

                // Wave speed
                const amrex::Real c = std::sqrt(grav_z * interior_h);

                // Calculation of Flather formula. "val" = velocity * h
                const auto Flather_val =
                    boundary_val + (ori.isLow() ? -1.0_rt : 1.0_rt) * c *
                                       (interior_h - boundary_h);

                // Use external (prescribed) velocity if inflow:
                // - Assesses inflow by the whole column, not the local value
                // - Assumes values are up-to-date from another fillpatch op
                const bool prescribed_inflow =
                    ori.isLow() ? boundary_val > 0.0_rt : boundary_val < 0.0_rt;

                // Clip the velocity in the interior to prevent inflow at
                // outflow; this is consistent with the line integral
                // calculations
                const auto local_internal_vel =
                    ori.isLow() ? amrex::min(ref_arr(iv_adj_cc, idir), 0.0_rt)
                                : amrex::max(ref_arr(iv_adj_cc, idir), 0.0_rt);

                // Normal outflow case:
                // *) For a given column, apply scale to the interior velocity,
                //    but only to fully liquid cells; leave mixed cells
                //    unchanged.

                // Edge cases:
                // 1) If the interior column contains no fully liquid cells,
                //    there is nothing to scale. Instead of scaling the internal
                //    profile, override the interior velocity with scaled
                //    external profile.
                // 2) During initialization, the scaling can be very aggressive.
                //    Plus, when the internal and external profiles are very
                //    different, the scaling can lead to rapid acceleration.
                //    Switch to the external profile when the changes are rapid.
                // 3) If the boundary velocity is 0, then the scaling based on
                //    external quantities is undefined. Use a regular outflow
                //    (Neumann) condition, even if edge cases 1) or 2) apply.
                // 4) If the boundary contains no liquid, the Flather velocity
                //    can be very large. Use a regular outflow (Neumann)
                //    condition.

                bool interior_valid =
                    std::abs(interior_liq) > v_threshold * interior_h;
                bool boundary_valid =
                    std::abs(boundary_val) > v_threshold * boundary_h;

                // Use boundary data by default (edge case 1)
                bool override_interior = true;
                auto scale_interior = 1.0_rt;
                if (interior_valid) {
                    // Normal case: calculate velocity scale with Flather
                    scale_interior =
                        (Flather_val - interior_mix) / interior_liq;

                    const bool interior_unbounded =
                        scale_interior > vscale_max ||
                        scale_interior < vscale_min;
                    // If scale is unbounded, override interior data
                    // (edge case 2)
                    override_interior = interior_unbounded;
                    // If scale is unbounded, revert to 1
                    // (part of edge case 3)
                    scale_interior =
                        interior_unbounded ? 1.0_rt : scale_interior;
                }

                // Never override interior with invalid boundary data
                // (part of edge case 3)
                override_interior = override_interior && boundary_valid;

                auto local_vel = 0.0_rt;
                auto scaled_vel = 0.0_rt;
                // Apply scale to velocity, depends on direction
                if (prescribed_inflow || override_interior) {
                    local_vel = arr(iv, fcomp);
                    scaled_vel = local_vel * (Flather_val / boundary_val);
                } else {
                    local_vel = local_internal_vel;
                    scaled_vel = local_vel * scale_interior;
                }

                // Only use if advecting liquid (or zero velocity):
                // The averages used to calculate the Flather velocity are only
                // performed on cells containing liquid; therefore, the scaling
                // of the velocity is only valid on those cells.
                const bool outflow =
                    ori.isLow() ? scaled_vel <= 0.0_rt : scaled_vel >= 0.0_rt;
                const bool inflow_any_liq = !outflow && boundary_vof > tiny;
                const bool outflow_only_liq =
                    outflow && interior_vof >= 1.0_rt - tiny;
                if (boundary_h > tiny && (outflow_only_liq || inflow_any_liq ||
                                          (outflow && override_interior))) {
                    arr(iv, fcomp) = scaled_vel;
                } else if (outflow) {
                    // edge case 4 (boundary_h <= tiny)
                    arr(iv, fcomp) = local_vel;
                }
            });
        }
    }
}

} // namespace kynema_sgf
