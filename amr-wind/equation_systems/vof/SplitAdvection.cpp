#include "amr-wind/equation_systems/vof/SplitAdvection.H"
#include "amr-wind/equation_systems/vof/split_advection.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

void multiphase::split_advection_step(
    int isweep,
    int iorder,
    int nlevels,
    Field& dof_field,
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> const& fluxes,
    ScratchField& fluxC,
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> const& advas,
    Field const& u_mac,
    Field const& v_mac,
    Field const& w_mac,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> BCs,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real time,
    amrex::Real dt,
    bool rm_debris)
{
    BL_PROFILE("amr-wind::multiphase::split_advection_step");

    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MFItInfo mfi_info;
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                .SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
             ++mfi) {
            const auto& bx = mfi.tilebox();
            amrex::FArrayBox tmpfab(
                amrex::grow(bx, 1), 2, amrex::The_Async_Arena());
            tmpfab.setVal<amrex::RunOn::Device>(0.0_rt);

            // Compression term coefficient
            if (iorder == 0) {
                multiphase::cmask_loop(
                    bx, dof_field(lev).array(mfi), fluxC(lev).array(mfi));
            }

            // Calculate fluxes involved in this stage of split advection
            multiphase::split_compute_fluxes(
                lev, bx, isweep + iorder, dof_field(lev).const_array(mfi),
                u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                w_mac(lev).const_array(mfi), (*advas[lev][0]).array(mfi),
                (*advas[lev][1]).array(mfi), (*advas[lev][2]).array(mfi),
                (*fluxes[lev][0]).array(mfi), (*fluxes[lev][1]).array(mfi),
                (*fluxes[lev][2]).array(mfi), BCs, tmpfab.dataPtr(), geom, dt);
        }
    }

    // Average down fluxes for current component
    const int dir = 2 - ((isweep + iorder) % 3);
    for (int lev = nlevels - 1; lev > 0; --lev) {
        amrex::IntVect rr =
            geom[lev].Domain().size() / geom[lev - 1].Domain().size();
        amrex::average_down_faces(
            *(fluxes[lev][dir]), (*fluxes[lev - 1][dir]), rr, geom[lev - 1]);
    }

    // Compute advection from first sweep
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MFItInfo mfi_info;
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                .SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
             ++mfi) {
            const auto& bx = mfi.tilebox();
            // Sum fluxes from this stage of advection
            multiphase::split_compute_sum(
                lev, bx, isweep + iorder, dof_field(lev).array(mfi),
                fluxC(lev).const_array(mfi), u_mac(lev).const_array(mfi),
                v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                (*fluxes[lev][0]).const_array(mfi),
                (*fluxes[lev][1]).const_array(mfi),
                (*fluxes[lev][2]).const_array(mfi), geom, dt);

            // Remove debris if desired after last step
            if (rm_debris && iorder == 2) {
                multiphase::debris_loop(bx, dof_field(lev).array(mfi));
            }
        }
    }

    // Average down vof and communicate
    dof_field.fillpatch(time);
}

void multiphase::split_compute_fluxes(
    const int lev,
    amrex::Box const& bx,
    const int isweep,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::Array4<amrex::Real> const& aax,
    amrex::Array4<amrex::Real> const& aay,
    amrex::Array4<amrex::Real> const& aaz,
    amrex::Array4<amrex::Real> const& fx,
    amrex::Array4<amrex::Real> const& fy,
    amrex::Array4<amrex::Real> const& fz,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> BCs,
    amrex::Real* p,
    amrex::Vector<amrex::Geometry> geom,
    const amrex::Real dt)
{
    BL_PROFILE("amr-wind::multiphase::split_compute_fluxes");
    amrex::Box const& bxg1 = amrex::grow(bx, 1);

    const amrex::Real dx = geom[lev].CellSize(0);
    const amrex::Real dy = geom[lev].CellSize(1);
    const amrex::Real dz = geom[lev].CellSize(2);
    amrex::Real dtdx = dt / dx;
    amrex::Real dtdy = dt / dy;
    amrex::Real dtdz = dt / dz;

    amrex::Box const& domain = geom[lev].Domain();
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    amrex::Array4<amrex::Real> vofL = makeArray4(p, bxg1, 1);
    p += vofL.size();
    amrex::Array4<amrex::Real> vofR = makeArray4(p, bxg1, 1);

    if (isweep % 3 == 0) {
        sweep_fluxes(2, bx, dtdz, wmac, volfrac, vofL, vofR);
        amrex::Box const& zbx = amrex::surroundingNodes(bx, 2);
        amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            fluxes_bc_save(
                i, j, k, 2, dt * wmac(i, j, k), fz, vofL, vofR, aaz, BCs,
                domlo.z, domhi.z);
        });
    } else if (isweep % 3 == 1) {
        sweep_fluxes(1, bx, dtdy, vmac, volfrac, vofL, vofR);
        amrex::Box const& ybx = amrex::surroundingNodes(bx, 1);
        amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            fluxes_bc_save(
                i, j, k, 1, dt * vmac(i, j, k), fy, vofL, vofR, aay, BCs,
                domlo.y, domhi.y);
        });
    } else {
        sweep_fluxes(0, bx, dtdx, umac, volfrac, vofL, vofR);
        amrex::Box const& xbx = amrex::surroundingNodes(bx, 0);
        amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            fluxes_bc_save(
                i, j, k, 0, dt * umac(i, j, k), fx, vofL, vofR, aax, BCs,
                domlo.x, domhi.x);
        });
    }
}

void multiphase::split_compute_sum(
    const int lev,
    amrex::Box const& bx,
    const int isweep,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxC,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    amrex::Vector<amrex::Geometry> geom,
    const amrex::Real dt)
{
    BL_PROFILE("amr-wind::multiphase::split_compute_sum");

    const auto dxinv = geom[lev].InvCellSizeArray();

    if (isweep % 3 == 0) {
        sweep_balance(2, dt, dxinv[2], bx, wmac, volfrac, fz, fluxC);
    } else if (isweep % 3 == 1) {
        sweep_balance(1, dt, dxinv[1], bx, vmac, volfrac, fy, fluxC);
    } else {
        sweep_balance(0, dt, dxinv[0], bx, umac, volfrac, fx, fluxC);
    }
}

void multiphase::cmask_loop(
    amrex::Box const& bx,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real> const& fluxC)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        multiphase::c_mask(i, j, k, volfrac, fluxC);
    });
}

void multiphase::debris_loop(
    amrex::Box const& bx, amrex::Array4<amrex::Real> const& volfrac)
{
    // Remove VOF debris
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        multiphase::remove_vof_debris(i, j, k, volfrac);
    });
}

void multiphase::sweep_fluxes(
    const int dir,
    amrex::Box const& bx,
    const amrex::Real dtdx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real> const& vofL,
    amrex::Array4<amrex::Real> const& vofR)
{

    BL_PROFILE("amr-wind::multiphase::sweep_fluxes");

    amrex::Box const& bxg1 = amrex::grow(bx, 1);

    int ii = (dir == 0) ? 1 : 0;
    int jj = (dir == 1) ? 1 : 0;
    int kk = (dir == 2) ? 1 : 0;

    amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::Real velL = vel_mac(i, j, k);
        amrex::Real velR = vel_mac(i + ii, j + jj, k + kk);
        multiphase::eulerian_implicit(
            i, j, k, dir, dtdx, velL, velR, volfrac, vofL, vofR);
    });
}

void multiphase::sweep_balance(
    const int dir,
    const amrex::Real dt,
    const amrex::Real dxi,
    amrex::Box const& bx,
    amrex::Array4<amrex::Real const> const& vel_mac,
    amrex::Array4<amrex::Real> const& volfrac,
    amrex::Array4<amrex::Real const> const& fluxF,
    amrex::Array4<amrex::Real const> const& fluxC)
{

    BL_PROFILE("amr-wind::multiphase::sweep_balance");

    int ii = (dir == 0) ? 1 : 0;
    int jj = (dir == 1) ? 1 : 0;
    int kk = (dir == 2) ? 1 : 0;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::Real velL = vel_mac(i, j, k) * dt * dxi;
        amrex::Real velR = vel_mac(i + ii, j + jj, k + kk) * dt * dxi;
        multiphase::balance_eulerian_fluxes(
            i, j, k, dir, dxi, velL, velR, volfrac, fluxF, fluxC);
    });
}

} // namespace amr_wind
