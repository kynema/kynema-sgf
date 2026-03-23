#ifndef VOF_MOMENTUM_FLUX_H
#define VOF_MOMENTUM_FLUX_H

#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::multiphase {

static AMREX_INLINE void hybrid_fluxes(
    const FieldRepo& repo,
    const int ncomp,
    const amrex::Gpu::DeviceVector<int>& iconserv,
    ScratchField& flux_x,
    ScratchField& flux_y,
    ScratchField& flux_z,
    const Field& dof_field,
    const Field& dof_nph,
    const Field& src_term,
    const Field& rho_o,
    const Field& rho_nph,
    const Field& u_mac,
    const Field& v_mac,
    const Field& w_mac,
    amrex::Vector<amrex::BCRec> const& velbc,
    amrex::BCRec const* velbc_d,
    amrex::Vector<amrex::BCRec> const& rhobc,
    amrex::BCRec const* rhobc_d,
    const amrex::Real dt,
    godunov::scheme mflux_scheme,
    bool allow_inflow_on_outflow,
    bool use_forces_in_trans,
    bool pre_multiplied_src_term = false)
{
    // Get geometry
    const auto& geom = repo.mesh().Geom();
    // Get advected alpha fields
    // At this point in the solve, they have been converted to advected density
    const auto& advrho_x = repo.get_field("advalpha_x");
    const auto& advrho_y = repo.get_field("advalpha_y");
    const auto& advrho_z = repo.get_field("advalpha_z");
    // Get VOF
    const auto& vof = repo.get_field("vof").state(amr_wind::FieldState::Old);

    // Create scratch arrays for local flux storage
    auto ftmp_x =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::XFACE);
    auto ftmp_y =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::YFACE);
    auto ftmp_z =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::ZFACE);

    auto face_x =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::XFACE);
    auto face_y =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::YFACE);
    auto face_z =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::ZFACE);

    // Create iconserv = 0 array for density, avoid multiplication by face vel
    amrex::Gpu::DeviceVector<int> idnsty;
    idnsty.resize(1, -1);

    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        // Set up temporary arrays for fluxes that will replace other quantities
        // form multifab for transport variable and source term
        amrex::MultiFab q(
            dof_field(lev).boxArray(), dof_field(lev).DistributionMap(), ncomp,
            fvm::Godunov::nghost_state);
        amrex::MultiFab::Copy(
            q, dof_field(lev), 0, 0, ncomp, fvm::Godunov::nghost_state);
        amrex::MultiFab fq(
            src_term(lev).boxArray(), src_term(lev).DistributionMap(), ncomp,
            fvm::Godunov::nghost_src);
        amrex::MultiFab::Copy(
            fq, src_term(lev), 0, 0, ncomp, fvm::Godunov::nghost_src);

        for (int idim = 0; idim < dof_field.num_comp(); ++idim) {
            amrex::MultiFab::Multiply(
                q, rho_o(lev), 0, idim, 1, fvm::Godunov::nghost_state);
            if (!pre_multiplied_src_term) {
                amrex::MultiFab::Multiply(
                    fq, rho_o(lev), 0, idim, 1, fvm::Godunov::nghost_src);
            }
        }

        amrex::MultiFab frho(
            src_term(lev).boxArray(), src_term(lev).DistributionMap(), 1,
            fvm::Godunov::nghost_src);
        frho.setVal(0.0_rt);

        // Set up temporary arrays for NPH at boundaries
        amrex::MultiFab q_nph(
            dof_field(lev).boxArray(), dof_field(lev).DistributionMap(), ncomp,
            fvm::Godunov::nghost_state);
        amrex::MultiFab::Copy(
            q_nph, dof_nph(lev), 0, 0, ncomp, fvm::Godunov::nghost_state);
        // Density at NPH is fully calculated a priori
        // Multiply to get momentum at NPH
        for (int idim = 0; idim < dof_field.num_comp(); ++idim) {
            amrex::MultiFab::Multiply(q_nph, rho_nph(lev), 0, idim, 1, 1);
        }

        std::string advection_type = "Godunov";
        int limiter_type = PPM::UPWIND; // default
        if (mflux_scheme == godunov::scheme::MINMOD) {
            limiter_type = PPM::MINMOD;
        } else if (mflux_scheme != godunov::scheme::UPWIND) {
            amrex::Abort("Unknown limiter type in vof_momentum_flux.H");
        }

        const bool known_edge_state = false;
        const bool godunov_use_ppm = true;
        const bool is_velocity = false;
        const bool fluxes_are_area_weighted = false;

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

            // Arrays for easy reference (w for working array)
            auto F_x = flux_x(lev).array(mfi);
            auto F_y = flux_y(lev).array(mfi);
            auto F_z = flux_z(lev).array(mfi);
            auto Fw_x = (*ftmp_x)(lev).array(mfi);
            auto Fw_y = (*ftmp_y)(lev).array(mfi);
            auto Fw_z = (*ftmp_z)(lev).array(mfi);

            // Temporary divu
            amrex::FArrayBox divufab(
                amrex::grow(bx, 1), 1, amrex::The_Async_Arena());
            divufab.setVal<amrex::RunOn::Device>(0.0_rt);
            const auto& divu = divufab.array();

            // Compute momentum flux quantities (multiplied by face velocity)
            HydroUtils::ComputeFluxesOnBoxFromState(
                bx, ncomp, mfi, q.const_array(mfi), q_nph.const_array(mfi),
                Fw_x, Fw_y, Fw_z, (*face_x)(lev).array(mfi),
                (*face_y)(lev).array(mfi), (*face_z)(lev).array(mfi),
                known_edge_state, u_mac(lev).const_array(mfi),
                v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi), divu,
                fq.const_array(mfi), geom[lev], dt, velbc, velbc_d,
                iconserv.data(), godunov_use_ppm, use_forces_in_trans,
                is_velocity, fluxes_are_area_weighted, advection_type,
                limiter_type, allow_inflow_on_outflow);

            const auto& bxg1 = amrex::grow(bx, 1);
            const auto& xbx = amrex::surroundingNodes(bx, 0);
            const auto& ybx = amrex::surroundingNodes(bx, 1);
            const auto& zbx = amrex::surroundingNodes(bx, 2);

            // Where interface is present, replace current momentum flux
            // quantities with those from mflux_scheme
            auto volfrac = vof(lev).const_array(mfi);
            amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                bool vf_i = multiphase::interface_band(i, j, k, volfrac);
                // Fluxes in each direction, check with flag too
                if (xbx.contains(i, j, k)) {
                    bool vf_nb =
                        multiphase::interface_band(i - 1, j, k, volfrac);
                    if (vf_i || vf_nb) {
                        for (int n = 0; n < ncomp; ++n) {
                            F_x(i, j, k, n) = Fw_x(i, j, k, n);
                        }
                    }
                }
                if (ybx.contains(i, j, k)) {
                    bool vf_nb =
                        multiphase::interface_band(i, j - 1, k, volfrac);
                    if (vf_i || vf_nb) {
                        for (int n = 0; n < ncomp; ++n) {
                            F_y(i, j, k, n) = Fw_y(i, j, k, n);
                        }
                    }
                }
                if (zbx.contains(i, j, k)) {
                    bool vf_nb =
                        multiphase::interface_band(i, j, k - 1, volfrac);
                    if (vf_i || vf_nb) {
                        for (int n = 0; n < ncomp; ++n) {
                            F_z(i, j, k, n) = Fw_z(i, j, k, n);
                        }
                    }
                }
            });

            // Compute density flux quantities (not multiplied by face velocity)
            HydroUtils::ComputeFluxesOnBoxFromState(
                bx, 1, mfi, rho_o(lev).const_array(mfi), Fw_x, Fw_y, Fw_z,
                (*face_x)(lev).array(mfi), (*face_y)(lev).array(mfi),
                (*face_z)(lev).array(mfi), known_edge_state,
                u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                w_mac(lev).const_array(mfi), divu, frho.const_array(mfi),
                geom[lev], dt, rhobc, rhobc_d, idnsty.data(), godunov_use_ppm,
                use_forces_in_trans, is_velocity, fluxes_are_area_weighted,
                advection_type, limiter_type, allow_inflow_on_outflow);

            // At this point in the solve, these are advected
            auto ar_x = advrho_x(lev).const_array(mfi);
            auto ar_y = advrho_y(lev).const_array(mfi);
            auto ar_z = advrho_z(lev).const_array(mfi);

            // When interface is present, divide by interpolated density to get
            // Favre-averaged flux value. Multiply all fluxes with advected
            // face density to create mass-consistent fluxes
            amrex::ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                bool vf_i = multiphase::interface_band(i, j, k, volfrac);
                // Fluxes in each direction, check with flag too
                if (xbx.contains(i, j, k)) {
                    bool vf_nb =
                        multiphase::interface_band(i - 1, j, k, volfrac);
                    for (int n = 0; n < ncomp; ++n) {
                        if (vf_i || vf_nb) {
                            F_x(i, j, k, n) /= Fw_x(i, j, k, 0);
                        }
                        F_x(i, j, k, n) *= ar_x(i, j, k);
                    }
                }
                if (ybx.contains(i, j, k)) {
                    bool vf_nb =
                        multiphase::interface_band(i, j - 1, k, volfrac);
                    for (int n = 0; n < ncomp; ++n) {
                        if (vf_i || vf_nb) {
                            F_y(i, j, k, n) /= Fw_y(i, j, k, 0);
                        }
                        F_y(i, j, k, n) *= ar_y(i, j, k);
                    }
                }
                if (zbx.contains(i, j, k)) {
                    bool vf_nb =
                        multiphase::interface_band(i, j, k - 1, volfrac);
                    for (int n = 0; n < ncomp; ++n) {
                        if (vf_i || vf_nb) {
                            F_z(i, j, k, n) /= Fw_z(i, j, k, 0);
                        }
                        F_z(i, j, k, n) *= ar_z(i, j, k);
                    }
                }
            });
        }
    }
}
} // namespace amr_wind::multiphase

#endif
