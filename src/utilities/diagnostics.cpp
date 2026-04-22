#include "amr-wind/incflo.H"
#include "diagnostics.H"
#include "constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

void amr_wind::diagnostics::make_mask_addend(
    amrex::MultiFab& mfab,
    const amrex::MultiFab& mfab_mask,
    const amrex::Real mask_val,
    const amrex::Real set_val)
{
    const auto& arr = mfab.arrays();
    const auto& arr_mask = mfab_mask.const_arrays();
    amrex::ParallelFor(
        mfab, mfab.n_grow, mfab.n_comp,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            arr[nbx](i, j, k, n) = std::abs(arr_mask[nbx](i, j, k) - mask_val) <
                                           constants::TIGHT_TOL
                                       ? 0.0_rt
                                       : set_val;
        });
}

void amr_wind::diagnostics::get_field_extrema(
    amrex::Real& field_max_val,
    amrex::Real& field_min_val,
    const amr_wind::Field& field,
    const int comp,
    const int ncomp,
    const int nghost)
{
    const int finest_level = field.repo().num_active_levels() - 1;

    amrex::Real max_val_lev{constants::LOW_NUM};
    amrex::Real min_val_lev{constants::LARGE_NUM};
    field_max_val = max_val_lev;
    field_min_val = min_val_lev;
    for (int lev = 0; lev <= finest_level; lev++) {
        for (int n = comp; n < comp + ncomp; ++n) {
            max_val_lev = field(lev).max(n, nghost);
            min_val_lev = field(lev).min(n, nghost);
            field_min_val = amrex::min<amrex::Real>(field_min_val, min_val_lev);
            field_max_val = amrex::max<amrex::Real>(field_max_val, max_val_lev);
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(field_max_val);
    amrex::ParallelDescriptor::ReduceRealMin(field_min_val);
}

bool amr_wind::diagnostics::get_field_extrema(
    amrex::Real& field_max_val,
    amrex::Real& field_min_val,
    const amr_wind::Field& field,
    const amr_wind::Field& field_mask,
    const amrex::Real mask_val,
    const int comp,
    const int ncomp,
    const int nghost)

{
    const int finest_level = field.repo().num_active_levels() - 1;

    amrex::Real max_val_lev{constants::LOW_NUM};
    amrex::Real min_val_lev{constants::LARGE_NUM};
    amrex::Real min_mask_lev = min_val_lev;
    amrex::Real mask_min = min_mask_lev;
    field_max_val = max_val_lev;
    field_min_val = min_val_lev;
    for (int lev = 0; lev <= finest_level; lev++) {
        amrex::MultiFab mask_addend_max(
            field_mask(lev).boxArray(), field_mask(lev).DistributionMap(),
            ncomp, nghost);
        amrex::MultiFab mask_addend_min(
            field_mask(lev).boxArray(), field_mask(lev).DistributionMap(),
            ncomp, nghost);
        amr_wind::diagnostics::make_mask_addend(
            mask_addend_max, field_mask(lev), mask_val, constants::LOW_NUM);
        amr_wind::diagnostics::make_mask_addend(
            mask_addend_min, field_mask(lev), mask_val, constants::LARGE_NUM);
        min_mask_lev = mask_addend_min.min(0, nghost);
        mask_min = amrex::min(mask_min, min_mask_lev);
        amrex::MultiFab::Add(
            mask_addend_max, field(lev), comp, 0, ncomp, nghost);
        amrex::MultiFab::Add(
            mask_addend_min, field(lev), comp, 0, ncomp, nghost);
        for (int n = 0; n < ncomp; ++n) {
            max_val_lev = mask_addend_max.max(n, nghost);
            min_val_lev = mask_addend_min.min(n, nghost);
            field_min_val = amrex::min<amrex::Real>(field_min_val, min_val_lev);
            field_max_val = amrex::max<amrex::Real>(field_max_val, max_val_lev);
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(field_max_val);
    amrex::ParallelDescriptor::ReduceRealMin(field_min_val);
    amrex::ParallelDescriptor::ReduceRealMax(mask_min);

    const bool found_mask_val = (mask_min < constants::TIGHT_TOL);
    return found_mask_val;
}

amrex::Real amr_wind::diagnostics::get_vel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir,
    const amrex::Real factor)
{
    return amrex::ReduceMax(
        vel, level_mask, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& vel_arr,
            amrex::Array4<int const> const& mask_arr) -> amrex::Real {
            amrex::Real max_fab = -1.0e8_rt;
            amrex::Loop(bx, [=, &max_fab](int i, int j, int k) {
                max_fab = amrex::max<amrex::Real>(
                    max_fab, mask_arr(i, j, k) > 0
                                 ? factor * vel_arr(i, j, k, vdir)
                                 : std::numeric_limits<amrex::Real>::lowest());
            });
            return max_fab;
        });
}

amrex::Real amr_wind::diagnostics::get_vel_max(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir)
{
    return get_vel_max(vel, level_mask, vdir, 1.0_rt);
}

amrex::Real amr_wind::diagnostics::get_vel_min(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir)
{
    return get_vel_max(vel, level_mask, vdir, -1.0_rt);
}

amrex::Real amr_wind::diagnostics::get_vel_loc(
    const amrex::MultiFab& vel,
    const amrex::iMultiFab& level_mask,
    const int vdir,
    const int ldir,
    amrex::Real vel_max,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    return amrex::ReduceMax(
        vel, level_mask, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& vel_arr,
            amrex::Array4<int const> const& mask_arr) -> amrex::Real {
            amrex::Real loc_fab = problo[ldir];
            amrex::Loop(bx, [=, &loc_fab](int i, int j, int k) {
                int idx = (ldir == 0 ? i : (ldir == 1 ? j : k));
                amrex::Real offset = 0.5_rt;
                amrex::Real loc = problo[ldir] + ((idx + offset) * dx[ldir]);
                bool mask_check = (mask_arr(i, j, k) > 0);
                bool loc_check =
                    (amrex::Math::abs(vel_max - vel_arr(i, j, k, vdir)) <
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
                loc_fab = amrex::max<amrex::Real>(
                    loc_fab, (mask_check && loc_check ? loc : problo[ldir]));
            });
            return loc_fab;
        });
}

amrex::Real amr_wind::diagnostics::get_macvel_max(
    const amrex::MultiFab& macvel,
    const amrex::iMultiFab& level_mask,
    const int vdir,
    const amrex::Real factor)
{
    return amrex::ReduceMax(
        macvel, level_mask, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& mvel_arr,
            amrex::Array4<int const> const& mask_arr) -> amrex::Real {
            amrex::Real max_fab = -1.0e8_rt;
            amrex::Loop(bx, [=, &max_fab](int i, int j, int k) {
                int ii = i - (vdir == 0 ? 1 : 0);
                int jj = j - (vdir == 1 ? 1 : 0);
                int kk = k - (vdir == 2 ? 1 : 0);
                max_fab = amrex::max<amrex::Real>(
                    max_fab, (mask_arr(i, j, k) + mask_arr(ii, jj, kk)) > 0
                                 ? factor * mvel_arr(i, j, k)
                                 : std::numeric_limits<amrex::Real>::lowest());
            });
            return max_fab;
        });
}

amrex::Real amr_wind::diagnostics::get_macvel_max(
    const amrex::MultiFab& macvel,
    const amrex::iMultiFab& level_mask,
    const int vdir)
{
    return get_macvel_max(macvel, level_mask, vdir, 1.0_rt);
}

amrex::Real amr_wind::diagnostics::get_macvel_min(
    const amrex::MultiFab& macvel,
    const amrex::iMultiFab& level_mask,
    const int vdir)
{
    return get_macvel_max(macvel, level_mask, vdir, -1.0_rt);
}

amrex::Real amr_wind::diagnostics::get_macvel_loc(
    const amrex::MultiFab& macvel,
    const amrex::iMultiFab& level_mask,
    const int vdir,
    const int ldir,
    amrex::Real mvel_max,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{

    return amrex::ReduceMax(
        macvel, level_mask, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& mvel_arr,
            amrex::Array4<int const> const& mask_arr) -> amrex::Real {
            amrex::Real loc_fab = problo[ldir];
            amrex::Loop(bx, [=, &loc_fab](int i, int j, int k) {
                int ii = i - (vdir == 0 ? 1 : 0);
                int jj = j - (vdir == 1 ? 1 : 0);
                int kk = k - (vdir == 2 ? 1 : 0);
                int idx = (ldir == 0 ? i : (ldir == 1 ? j : k));
                amrex::Real offset = (ldir == vdir ? 0.0_rt : 0.5_rt);
                amrex::Real loc = problo[ldir] + ((idx + offset) * dx[ldir]);
                bool mask_check =
                    (mask_arr(i, j, k) + mask_arr(ii, jj, kk) > 0);
                bool loc_check =
                    (amrex::Math::abs(mvel_max - mvel_arr(i, j, k)) <
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
                loc_fab = amrex::max<amrex::Real>(
                    loc_fab, (mask_check && loc_check ? loc : problo[ldir]));
            });
            return loc_fab;
        });
}

amrex::Array<amrex::Real, 24> amr_wind::diagnostics::PrintMaxVelLocations(
    const amr_wind::FieldRepo& repo, const std::string& header)
{
    BL_PROFILE("amr-wind::diagnostics::PrintMaxVelLocations");

    // Get fields
    const auto& vel = repo.get_field("velocity");
    const int finest_level = repo.num_active_levels() - 1;

    // Get infinity norm of velocities
    amrex::Real u_max{-1.0e8_rt}, v_max{-1.0e8_rt}, w_max{-1.0e8_rt};
    // Minima will be negated later
    amrex::Real u_min{-1.0e8_rt}, v_min{-1.0e8_rt}, w_min{-1.0e8_rt};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                repo.mesh().boxArray(lev + 1), repo.mesh().refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev), 1,
                0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        u_max = amrex::max<amrex::Real>(
            u_max, amr_wind::diagnostics::get_vel_max(vel(lev), level_mask, 0));

        u_min = amrex::max<amrex::Real>(
            u_min, amr_wind::diagnostics::get_vel_min(vel(lev), level_mask, 0));

        v_max = amrex::max<amrex::Real>(
            v_max, amr_wind::diagnostics::get_vel_max(vel(lev), level_mask, 1));

        v_min = amrex::max<amrex::Real>(
            v_min, amr_wind::diagnostics::get_vel_min(vel(lev), level_mask, 1));

        w_max = amrex::max<amrex::Real>(
            w_max, amr_wind::diagnostics::get_vel_max(vel(lev), level_mask, 2));

        w_min = amrex::max<amrex::Real>(
            w_min, amr_wind::diagnostics::get_vel_min(vel(lev), level_mask, 2));
    }

    // Do additional parallelism stuff
    amrex::ParallelDescriptor::ReduceRealMax(u_max);
    amrex::ParallelDescriptor::ReduceRealMax(v_max);
    amrex::ParallelDescriptor::ReduceRealMax(w_max);
    amrex::ParallelDescriptor::ReduceRealMax(u_min);
    amrex::ParallelDescriptor::ReduceRealMax(v_min);
    amrex::ParallelDescriptor::ReduceRealMax(w_min);

    // Negate minima
    u_min *= -1.0_rt;
    v_min *= -1.0_rt;
    w_min *= -1.0_rt;

    // Get locations of these extrema
    auto problo = (repo.mesh().Geom())[0].ProbLoArray();
    amrex::GpuArray<amrex::Real, 3> u_max_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> v_max_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> w_max_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> u_min_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> v_min_loc{problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> w_min_loc{problo[0], problo[1], problo[2]};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                repo.mesh().boxArray(lev + 1), repo.mesh().refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev), 1,
                0, amrex::MFInfo());
            level_mask.setVal(1);
        }

        problo = (repo.mesh().Geom())[lev].ProbLoArray();
        auto dx = (repo.mesh().Geom())[lev].CellSizeArray();

        // Loop coordinates directions
        for (int n = 0; n < 3; n++) {
            u_max_loc[n] = amrex::max<amrex::Real>(
                u_max_loc[n],
                amr_wind::diagnostics::get_vel_loc(
                    vel(lev), level_mask, 0, n, u_max, problo, dx));

            u_min_loc[n] = amrex::max<amrex::Real>(
                u_min_loc[n],
                amr_wind::diagnostics::get_vel_loc(
                    vel(lev), level_mask, 0, n, u_min, problo, dx));

            v_max_loc[n] = amrex::max<amrex::Real>(
                v_max_loc[n],
                amr_wind::diagnostics::get_vel_loc(
                    vel(lev), level_mask, 1, n, v_max, problo, dx));

            v_min_loc[n] = amrex::max<amrex::Real>(
                v_min_loc[n],
                amr_wind::diagnostics::get_vel_loc(
                    vel(lev), level_mask, 1, n, v_min, problo, dx));

            w_max_loc[n] = amrex::max<amrex::Real>(
                w_max_loc[n],
                amr_wind::diagnostics::get_vel_loc(
                    vel(lev), level_mask, 2, n, w_max, problo, dx));

            w_min_loc[n] = amrex::max<amrex::Real>(
                w_min_loc[n],
                amr_wind::diagnostics::get_vel_loc(
                    vel(lev), level_mask, 2, n, w_min, problo, dx));
        }
    }

    // Additional parallelism
    for (int n = 0; n < 3; ++n) {
        amrex::ParallelDescriptor::ReduceRealMax(u_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(v_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(w_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(u_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(v_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(w_min_loc[n]);
    }

    // Output results
    amrex::Print() << "\nL-inf norm vels: " << header << '\n'
                   << "........................................................"
                      "......................"
                   << '\n';

    amrex::Print() << "Max u: " << std::setw(20) << std::right << u_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << u_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_max_loc[2] << '\n';
    amrex::Print() << "Min u: " << std::setw(20) << std::right << u_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << u_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << u_min_loc[2] << '\n';

    amrex::Print() << "Max v: " << std::setw(20) << std::right << v_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << v_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_max_loc[2] << '\n';
    amrex::Print() << "Min v: " << std::setw(20) << std::right << v_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << v_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << v_min_loc[2] << '\n';

    amrex::Print() << "Max w: " << std::setw(20) << std::right << w_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << w_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_max_loc[2] << '\n';
    amrex::Print() << "Min w: " << std::setw(20) << std::right << w_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << w_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << w_min_loc[2] << '\n';

    amrex::Print() << "........................................................"
                      "......................"
                   << '\n'
                   << '\n';

    // Return array of answers (for testing)
    return amrex::Array<amrex::Real, 24>{
        u_max,        u_max_loc[0], u_max_loc[1], u_max_loc[2], u_min,
        u_min_loc[0], u_min_loc[1], u_min_loc[2], v_max,        v_max_loc[0],
        v_max_loc[1], v_max_loc[2], v_min,        v_min_loc[0], v_min_loc[1],
        v_min_loc[2], w_max,        w_max_loc[0], w_max_loc[1], w_max_loc[2],
        w_min,        w_min_loc[0], w_min_loc[1], w_min_loc[2]};
}

amrex::Array<amrex::Real, 24> amr_wind::diagnostics::PrintMaxMACVelLocations(
    const amr_wind::FieldRepo& repo, const std::string& header)
{
    BL_PROFILE("amr-wind::diagnostics::PrintMaxMACVelLocations");

    // Get fields
    const auto& u_mac = repo.get_field("u_mac");
    const auto& v_mac = repo.get_field("v_mac");
    const auto& w_mac = repo.get_field("w_mac");
    const int finest_level = repo.num_active_levels() - 1;

    // Get infinity norm of mac velocities
    amrex::Real uMAC_max{-1.0e8_rt}, vMAC_max{-1.0e8_rt}, wMAC_max{-1.0e8_rt};
    // Minima will be negated later
    amrex::Real uMAC_min{-1.0e8_rt}, vMAC_min{-1.0e8_rt}, wMAC_min{-1.0e8_rt};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        // Do it with a ghost cell for the sake of checking faces
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            // MultiFab with ghost cell
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev), 1,
                1, amrex::MFInfo());
            // Set to 0 (ghost cells will be 0)
            level_mask.setVal(0);
            // Populate non-ghost cells with fine mask
            amrex::iMultiFab::Copy(
                level_mask,
                makeFineMask(
                    repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                    repo.mesh().boxArray(lev + 1), repo.mesh().refRatio(lev), 1,
                    0),
                0, 0, 1, amrex::IntVect(0));
        } else {
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev), 1,
                1, amrex::MFInfo());
            level_mask.setVal(1);
        }

        uMAC_max = amrex::max<amrex::Real>(
            uMAC_max,
            amr_wind::diagnostics::get_macvel_max(u_mac(lev), level_mask, 0));

        uMAC_min = amrex::max<amrex::Real>(
            uMAC_min,
            amr_wind::diagnostics::get_macvel_min(u_mac(lev), level_mask, 0));

        vMAC_max = amrex::max<amrex::Real>(
            vMAC_max,
            amr_wind::diagnostics::get_macvel_max(v_mac(lev), level_mask, 1));

        vMAC_min = amrex::max<amrex::Real>(
            vMAC_min,
            amr_wind::diagnostics::get_macvel_min(v_mac(lev), level_mask, 1));

        wMAC_max = amrex::max<amrex::Real>(
            wMAC_max,
            amr_wind::diagnostics::get_macvel_max(w_mac(lev), level_mask, 2));

        wMAC_min = amrex::max<amrex::Real>(
            wMAC_min,
            amr_wind::diagnostics::get_macvel_min(w_mac(lev), level_mask, 2));
    }

    // Do additional parallelism stuff
    amrex::ParallelDescriptor::ReduceRealMax(uMAC_max);
    amrex::ParallelDescriptor::ReduceRealMax(vMAC_max);
    amrex::ParallelDescriptor::ReduceRealMax(wMAC_max);
    amrex::ParallelDescriptor::ReduceRealMax(uMAC_min);
    amrex::ParallelDescriptor::ReduceRealMax(vMAC_min);
    amrex::ParallelDescriptor::ReduceRealMax(wMAC_min);

    // Negate minima
    uMAC_min *= -1.0_rt;
    vMAC_min *= -1.0_rt;
    wMAC_min *= -1.0_rt;

    // Get locations of these extrema
    auto problo = (repo.mesh().Geom())[0].ProbLoArray();
    amrex::GpuArray<amrex::Real, 3> uMAC_max_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> vMAC_max_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> wMAC_max_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> uMAC_min_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> vMAC_min_loc{
        problo[0], problo[1], problo[2]};
    amrex::GpuArray<amrex::Real, 3> wMAC_min_loc{
        problo[0], problo[1], problo[2]};
    for (int lev = 0; lev <= finest_level; lev++) {
        // Use level_mask to only count finest level present
        // Do it with a ghost cell for the sake of checking faces
        amrex::iMultiFab level_mask;
        if (lev < finest_level) {
            // MultiFab with ghost cell
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev), 1,
                1, amrex::MFInfo());
            // Set to 0 (ghost cells will be 0)
            level_mask.setVal(0);
            // Populate non-ghost cells with fine mask
            amrex::iMultiFab::Copy(
                level_mask,
                makeFineMask(
                    repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev),
                    repo.mesh().boxArray(lev + 1), repo.mesh().refRatio(lev), 1,
                    0),
                0, 0, 1, amrex::IntVect(0));
        } else {
            level_mask.define(
                repo.mesh().boxArray(lev), repo.mesh().DistributionMap(lev), 1,
                1, amrex::MFInfo());
            level_mask.setVal(1);
        }

        problo = (repo.mesh().Geom())[lev].ProbLoArray();
        auto dx = (repo.mesh().Geom())[lev].CellSizeArray();

        // Loop coordinates directions
        for (int n = 0; n < 3; n++) {
            uMAC_max_loc[n] = amrex::max<amrex::Real>(
                uMAC_max_loc[n],
                amr_wind::diagnostics::get_macvel_loc(
                    u_mac(lev), level_mask, 0, n, uMAC_max, problo, dx));

            uMAC_min_loc[n] = amrex::max<amrex::Real>(
                uMAC_min_loc[n],
                amr_wind::diagnostics::get_macvel_loc(
                    u_mac(lev), level_mask, 0, n, uMAC_min, problo, dx));

            vMAC_max_loc[n] = amrex::max<amrex::Real>(
                vMAC_max_loc[n],
                amr_wind::diagnostics::get_macvel_loc(
                    v_mac(lev), level_mask, 1, n, vMAC_max, problo, dx));

            vMAC_min_loc[n] = amrex::max<amrex::Real>(
                vMAC_min_loc[n],
                amr_wind::diagnostics::get_macvel_loc(
                    v_mac(lev), level_mask, 1, n, vMAC_min, problo, dx));

            wMAC_max_loc[n] = amrex::max<amrex::Real>(
                wMAC_max_loc[n],
                amr_wind::diagnostics::get_macvel_loc(
                    w_mac(lev), level_mask, 2, n, wMAC_max, problo, dx));

            wMAC_min_loc[n] = amrex::max<amrex::Real>(
                wMAC_min_loc[n],
                amr_wind::diagnostics::get_macvel_loc(
                    w_mac(lev), level_mask, 2, n, wMAC_min, problo, dx));
        }
    }

    // Additional parallelism
    for (int n = 0; n < 3; ++n) {
        amrex::ParallelDescriptor::ReduceRealMax(uMAC_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(vMAC_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(wMAC_max_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(uMAC_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(vMAC_min_loc[n]);
        amrex::ParallelDescriptor::ReduceRealMax(wMAC_min_loc[n]);
    }

    // Output results
    amrex::Print() << "\nL-inf norm MAC vels: " << header << '\n'
                   << "........................................................"
                      "......................"
                   << '\n';

    amrex::Print() << "Max u: " << std::setw(20) << std::right << uMAC_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << uMAC_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_max_loc[2] << '\n';
    amrex::Print() << "Min u: " << std::setw(20) << std::right << uMAC_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << uMAC_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << uMAC_min_loc[2] << '\n';

    amrex::Print() << "Max v: " << std::setw(20) << std::right << vMAC_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << vMAC_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_max_loc[2] << '\n';
    amrex::Print() << "Min v: " << std::setw(20) << std::right << vMAC_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << vMAC_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << vMAC_min_loc[2] << '\n';

    amrex::Print() << "Max w: " << std::setw(20) << std::right << wMAC_max;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << wMAC_max_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_max_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_max_loc[2] << '\n';
    amrex::Print() << "Min w: " << std::setw(20) << std::right << wMAC_min;
    amrex::Print() << " |  Location (x,y,z): ";
    amrex::Print() << std::setw(10) << std::right << wMAC_min_loc[0] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_min_loc[1] << ", ";
    amrex::Print() << std::setw(10) << std::right << wMAC_min_loc[2] << '\n';

    amrex::Print() << "........................................................"
                      "......................"
                   << '\n'
                   << '\n';

    // Return array of answers (for testing)
    return amrex::Array<amrex::Real, 24>{
        uMAC_max, uMAC_max_loc[0], uMAC_max_loc[1], uMAC_max_loc[2],
        uMAC_min, uMAC_min_loc[0], uMAC_min_loc[1], uMAC_min_loc[2],
        vMAC_max, vMAC_max_loc[0], vMAC_max_loc[1], vMAC_max_loc[2],
        vMAC_min, vMAC_min_loc[0], vMAC_min_loc[1], vMAC_min_loc[2],
        wMAC_max, wMAC_max_loc[0], wMAC_max_loc[1], wMAC_max_loc[2],
        wMAC_min, wMAC_min_loc[0], wMAC_min_loc[1], wMAC_min_loc[2]};
}

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(const std::string& header)
{
    BL_PROFILE("amr-wind::incflo::PrintMaxValues");

    amrex::Print() << "\nL-inf norm summary: " << header << '\n'
                   << "........................................................"
                      "......................";

    for (int lev = 0; lev <= finest_level; lev++) {
        amrex::Print() << "\nLevel " << lev << '\n';

        const auto& vel = icns().fields().field;
        amrex::Print() << "  " << std::setw(16) << std::left << vel.name();
        for (int i = 0; i < vel.num_comp(); ++i) {
            amrex::Print() << std::setw(20) << std::right << vel(lev).norm0(i);
        }
        amrex::Print() << '\n';

        const auto& gradp = grad_p();
        amrex::Print() << "  " << std::setw(16) << std::left << gradp.name();
        for (int i = 0; i < gradp.num_comp(); ++i) {
            amrex::Print() << std::setw(20) << std::right
                           << gradp(lev).norm0(i);
        }
        amrex::Print() << '\n';

        for (auto& eqn : scalar_eqns()) {
            auto& field = eqn->fields().field;
            amrex::Print() << "  " << std::setw(16) << std::left
                           << field.name();
            for (int i = 0; i < field.num_comp(); ++i) {
                amrex::Print()
                    << std::setw(20) << std::right << field(lev).norm0(i);
            }
            amrex::Print() << '\n';
        }
    }
    amrex::Print() << "........................................................"
                      "......................"
                   << '\n'
                   << '\n';
}

void incflo::PrintMaxVelLocations(const std::string& header)
{
    amr_wind::diagnostics::PrintMaxVelLocations(repo(), header);
}

//
// Print the maximum values of the velocity components and velocity divergence
//
void incflo::PrintMaxVel(int lev) const
{
    BL_PROFILE("amr-wind::incflo::PrintMaxVel");
    amrex::Print() << "max(abs(u/v/w))  = " << velocity()(lev).norm0(0) << "  "
                   << velocity()(lev).norm0(1) << "  "
                   << velocity()(lev).norm0(2) << "  " << '\n';
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev) const
{
    BL_PROFILE("amr-wind::incflo::PrintMaxGp");
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = " << grad_p()(lev).norm0(0)
                   << "  " << grad_p()(lev).norm0(1) << "  "
                   << grad_p()(lev).norm0(2) << "  " << pressure()(lev).norm0(0)
                   << "  " << '\n';
}

void incflo::CheckForNans(int lev) const
{
    BL_PROFILE("amr-wind::incflo::CheckForNans");
    bool ro_has_nans = density()(lev).contains_nan(false);
    bool ug_has_nans = velocity()(lev).contains_nan(false);
    bool vg_has_nans = velocity()(lev).contains_nan(true);
    bool wg_has_nans = velocity()(lev).contains_nan(true);
    bool pg_has_nans = pressure()(lev).contains_nan(false);

    if (ro_has_nans) {
        amrex::Print() << "WARNING: ro contains NaNs!!!";
    }

    if (ug_has_nans) {
        amrex::Print() << "WARNING: u contains NaNs!!!";
    }

    if (vg_has_nans) {
        amrex::Print() << "WARNING: v contains NaNs!!!";
    }

    if (wg_has_nans) {
        amrex::Print() << "WARNING: w contains NaNs!!!";
    }

    if (pg_has_nans) {
        amrex::Print() << "WARNING: p contains NaNs!!!";
    }
}
