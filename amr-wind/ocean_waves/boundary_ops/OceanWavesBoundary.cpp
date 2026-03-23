#include "amr-wind/CFDSim.H"
#include "amr-wind/ocean_waves/boundary_ops/OceanWavesBoundary.H"
#include "amr-wind/ocean_waves/boundary_ops/OceanWavesFillInflow.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/constants.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

OceanWavesBoundary::OceanWavesBoundary(CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_ow_velocity(sim.repo().get_field("ow_velocity"))
    , m_ow_vof(sim.repo().get_field("ow_vof"))
{
    // Check for if boundary planes are being used; disable if so
    if (sim.physics_manager().contains("ABL")) {
        if (sim.physics_manager().get<amr_wind::ABL>().bndry_plane().mode() ==
            io_mode::input) {
            // Turn off ow_bndry; will rely on bndry_plane for fills
            m_activate_ow_bndry = false;
        }
        if (sim.physics_manager().get<amr_wind::ABL>().abl_mpl().is_active()) {
            amrex::Abort(
                "OceanWavesBoundary: not currently compatible with ABL MPL "
                "implementation.");
        }
    }
    // Get liquid density, will only be used if vof is present
    if (sim.physics_manager().contains("MultiPhase")) {
        m_rho1 = sim.physics_manager().get<amr_wind::MultiPhase>().rho1();
    }
}

void OceanWavesBoundary::post_init_actions()
{
    BL_PROFILE("amr-wind::OceanWavesBoundary::post_init_actions");
    if (m_activate_ow_bndry) {
        m_repo.get_field("velocity")
            .register_fill_patch_op<OceanWavesFillInflow>(
                m_mesh, m_time, *this);
        m_vof_exists = m_repo.field_exists("vof");
        if (m_vof_exists) {
            m_repo.get_field("vof")
                .register_fill_patch_op<OceanWavesFillInflow>(
                    m_mesh, m_time, *this);
            m_repo.get_field("density")
                .register_fill_patch_op<OceanWavesFillInflow>(
                    m_mesh, m_time, *this);
        }

        m_terrain_exists = m_repo.int_field_exists("terrain_blank");
        if (m_terrain_exists) {
            m_terrain_blank_ptr = &m_repo.get_int_field("terrain_blank");
        }
    }
}

void OceanWavesBoundary::set_velocity(
    const int lev,
    const amrex::Real /*time*/,
    const Field& fld,
    amrex::MultiFab& mfab,
    const int dcomp,
    const int orig_comp) const
{

    if (!m_activate_ow_bndry) {
        return;
    }

    BL_PROFILE("amr-wind::OceanWavesBoundary::set_velocity");

    const auto& geom = m_mesh.Geom(lev);
    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    const bool terrain_and_vof_exist = m_terrain_exists && m_vof_exists;

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow) &&
            (bctype[ori] != BC::wave_generation)) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

        auto shift_to_interior =
            amrex::IntVect::TheDimensionVector(idir) * (ori.isLow() ? 1 : -1);

        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            amrex::IntVect shift_to_cc = {0, 0, 0};
            const auto& bx = utils::face_aware_boundary_box_intersection(
                shift_to_cc, gbx, dbx, ori);
            if (!bx.ok()) {
                continue;
            }

            const auto& targ_vof = m_ow_vof(lev).const_array(mfi);
            const auto& targ_arr = m_ow_velocity(lev).const_array(mfi);
            const auto& arr = mfab[mfi].array();
            const int numcomp = mfab.nComp();

            const auto terrain_blank_flags =
                terrain_and_vof_exist
                    ? (*m_terrain_blank_ptr)(lev).const_array(mfi)
                    : amrex::Array4<int const>();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::IntVect iv{i, j, k};
                for (int n = 0; n < numcomp; n++) {
                    if (targ_vof(iv) > constants::TIGHT_TOL) {
                        arr(iv, dcomp + n) = targ_arr(iv, orig_comp + n);
                    }
                    if (terrain_and_vof_exist) {
                        // Terrain-blanked adjacent cell means 0 velocity
                        if (terrain_blank_flags(
                                iv + shift_to_cc + shift_to_interior) == 1) {
                            arr(iv, dcomp + n) = 0.0_rt;
                        }
                    }
                }
            });
        }
    }
}

void OceanWavesBoundary::set_vof(
    const int lev,
    const amrex::Real /*time*/,
    const Field& fld,
    amrex::MultiFab& mfab) const
{

    if (!m_activate_ow_bndry) {
        return;
    }

    BL_PROFILE("amr-wind::OceanWavesBoundary::set_vof");

    const auto& geom = m_mesh.Geom(lev);
    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow) &&
            (bctype[ori] != BC::wave_generation)) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            const auto& bx =
                utils::face_aware_boundary_box_intersection(gbx, dbx, ori);
            if (!bx.ok()) {
                continue;
            }

            const auto& targ_arr = m_ow_vof(lev).const_array(mfi);
            const auto& arr = mfab[mfi].array();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                arr(i, j, k) = targ_arr(i, j, k);
            });
        }
    }
}

void OceanWavesBoundary::set_density(
    const int lev,
    const amrex::Real /*time*/,
    const Field& fld,
    amrex::MultiFab& mfab) const
{

    if (!m_activate_ow_bndry || m_rho1 < 0.0_rt) {
        return;
    }

    BL_PROFILE("amr-wind::OceanWavesBoundary::set_density");

    const auto& geom = m_mesh.Geom(lev);
    const auto& bctype = fld.bc_type();
    const int nghost = 1;
    const amrex::Real rho1 = m_rho1;
    const auto& domain = geom.growPeriodicDomain(nghost);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow) &&
            (bctype[ori] != BC::wave_generation)) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& dbx = ori.isLow() ? amrex::adjCellLo(domain, idir, nghost)
                                      : amrex::adjCellHi(domain, idir, nghost);

        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            auto gbx = amrex::grow(mfi.validbox(), nghost);
            const auto& bx =
                utils::face_aware_boundary_box_intersection(gbx, dbx, ori);
            if (!bx.ok()) {
                continue;
            }

            const auto& targ_vof = m_ow_vof(lev).const_array(mfi);
            const auto& arr = mfab[mfi].array();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                // Assume density is correct for gas phase only
                arr(i, j, k) = (targ_vof(i, j, k) * rho1) +
                               ((1.0_rt - targ_vof(i, j, k)) * arr(i, j, k));
            });
        }
    }
}

void OceanWavesBoundary::set_inflow_sibling_velocity(
    const int lev,
    const amrex::Real /*time*/,
    const Field& fld,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mfabs) const
{

    if (!m_activate_ow_bndry) {
        return;
    }

    BL_PROFILE("amr-wind::OceanWavesBoundary::set_inflow_sibling_velocity");

    const bool terrain_and_vof_exist = m_terrain_exists && m_vof_exists;
    const auto& bctype = fld.bc_type();
    const auto& geom = fld.repo().mesh().Geom(lev);

    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        if ((bctype[ori] != BC::mass_inflow) &&
            (bctype[ori] != BC::mass_inflow_outflow) &&
            (bctype[ori] != BC::wave_generation)) {
            continue;
        }

        const int idir = ori.coordDir();
        const auto& domain_box = geom.Domain();

        auto shift_to_interior =
            amrex::IntVect::TheDimensionVector(idir) * (ori.isLow() ? 1 : -1);

        for (int fdir = 0; fdir < AMREX_SPACEDIM; ++fdir) {

            // Only face-normal velocities populated here
            if (idir != fdir) {
                continue;
            }
            const auto& dbx = ori.isLow() ? amrex::bdryLo(domain_box, idir)
                                          : amrex::bdryHi(domain_box, idir);

            // Shift from valid face index to first cell-centered ghost
            amrex::IntVect shift_to_cc = {0, 0, 0};
            if (ori.isLow()) {
                --shift_to_cc[fdir];
            }

            auto& mfab = *mfabs[fdir];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                const auto& bx = vbx & dbx;
                if (!bx.ok()) {
                    continue;
                }

                const auto& targ_vof = m_ow_vof(lev).const_array(mfi);
                const auto& targ_arr = m_ow_velocity(lev).const_array(mfi);
                const auto& marr = mfab[mfi].array();

                const auto terrain_blank_flags =
                    terrain_and_vof_exist
                        ? (*m_terrain_blank_ptr)(lev).const_array(mfi)
                        : amrex::Array4<int const>();

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(
                            const int i, const int j, const int k) {
                        amrex::IntVect cc_iv = {i, j, k};
                        cc_iv += shift_to_cc;

                        if (targ_vof(cc_iv) > constants::TIGHT_TOL) {
                            marr(i, j, k, 0) = targ_arr(cc_iv, fdir);
                        }
                        if (terrain_and_vof_exist) {
                            // Terrain-blanked boundary-adjacent cell should set
                            // boundary velocity to 0
                            if (terrain_blank_flags(
                                    cc_iv + shift_to_interior) == 1) {
                                marr(i, j, k, 0) = 0.0_rt;
                            }
                        }
                    });
            }
        }
    }
}

} // namespace amr_wind
