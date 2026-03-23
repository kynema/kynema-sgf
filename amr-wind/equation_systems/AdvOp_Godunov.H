#ifndef ADVOP_GODUNOV_H
#define ADVOP_GODUNOV_H

#include <type_traits>

#include "amr-wind/convection/Godunov.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"

#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "hydro_utils.H"

#include "amr-wind/equation_systems/vof/vof_momentum_flux.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/** Godunov advection operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct AdvectionOp<
    PDE,
    fvm::Godunov,
    typename std::enable_if_t<std::is_base_of_v<ScalarTransport, PDE>>>
{
    AdvectionOp(
        CFDSim& /* sim */,
        PDEFields& fields_in,
        bool /* has_overset */,
        bool /* variable density */,
        bool /* mesh mapping */,
        bool /* is_anelastic */)
        : fields(fields_in)
        , density(fields_in.repo.get_field("density"))
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
    {
        amrex::ParmParse pp("incflo");
        pp.query("godunov_type", godunov_type);
        pp.query("godunov_use_forces_in_trans", godunov_use_forces_in_trans);
        if (pp.contains("use_ppm") || pp.contains("use_limiter")) {
            amrex::Abort(
                "Godunov: use_ppm and use_limiter are deprecated. Please "
                "update input file");
        }

        if (amrex::toLower(godunov_type) == "plm") {
            godunov_scheme = godunov::scheme::PLM;
        } else if (amrex::toLower(godunov_type) == "ppm") {
            godunov_scheme = godunov::scheme::PPM;
        } else if (amrex::toLower(godunov_type) == "ppm_nolim") {
            godunov_scheme = godunov::scheme::PPM_NOLIM;
            amrex::Print() << "WARNING: Using advection type ppm_nolim is not "
                              "recommended. Prefer using weno_z."
                           << '\n';
        } else if (amrex::toLower(godunov_type) == "bds") {
            godunov_scheme = godunov::scheme::BDS;
            advection_type = "BDS";
        } else if (
            amrex::toLower(godunov_type) == "weno" ||
            amrex::toLower(godunov_type) == "weno_js") {
            godunov_scheme = godunov::scheme::WENO_JS;
        } else if (amrex::toLower(godunov_type) == "weno_z") {
            godunov_scheme = godunov::scheme::WENOZ;
        } else {
            amrex::Abort(
                "Invalid godunov_type specified. For godunov_type select "
                "between plm, ppm, ppm_nolim, bds, weno_js, and weno_z. If no "
                "godunov_type is specified, the default weno_z is used.");
        }

        // Flux calculation used in multiphase portions of domain
        pp.query("mflux_type", mflux_type);
        if (amrex::toLower(mflux_type) == "minmod") {
            mflux_scheme = godunov::scheme::MINMOD;
        } else if (amrex::toLower(mflux_type) == "upwind") {
            mflux_scheme = godunov::scheme::UPWIND;
        } else {
            amrex::Abort("Invalid argument entered for mflux_type.");
        }

        amrex::ParmParse pp_eq{fields_in.field.base_name()};
        pp_eq.query(
            "allow_inflow_at_pressure_outflow", m_allow_inflow_on_outflow);

        // TODO: Need iconserv flag to be adjusted???
        iconserv.resize(PDE::ndim, 1);
    }

    void preadvect(
        const FieldState /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/)
    {}

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        static_assert(
            PDE::ndim == 1, "Invalid number of components for scalar");
        auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        // cppcheck-suppress constVariableReference
        auto& conv_term = fields.conv_term;
        const auto& dof_field = fields.field.state(fstate);
        const auto& dof_nph = fields.field.state(amr_wind::FieldState::NPH);

        auto flux_x =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto flux_y =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto flux_z =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::ZFACE);
        auto face_x =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto face_y =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto face_z =
            repo.create_scratch_field(PDE::ndim, 0, amr_wind::FieldLoc::ZFACE);

        // only needed if multiplying by rho below
        const auto& den = density.state(fstate);
        const auto& den_nph = density.state(amr_wind::FieldState::NPH);

        const bool mphase_vof = repo.field_exists("vof");

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
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
                auto rho_arr = den(lev).array(mfi);
                auto rho_nph_arr = den_nph(lev).array(mfi);
                auto tra_arr = dof_field(lev).array(mfi);
                auto tra_nph_arr = dof_nph(lev).array(mfi);
                amrex::FArrayBox rhotracfab;
                amrex::Array4<amrex::Real> rhotrac;
                amrex::FArrayBox rhotrac_nph_fab;
                amrex::Array4<amrex::Real> rhotrac_nph;
                amrex::FArrayBox srctrac_over_rho_fab;
                amrex::Array4<amrex::Real> srctrac_over_rho;
                const auto& src_arr = src_term(lev).const_array(mfi);

                if (PDE::multiply_rho) {
                    auto rhotrac_box =
                        amrex::grow(bx, fvm::Godunov::nghost_state);
                    auto srctrac_box =
                        amrex::grow(bx, fvm::Godunov::nghost_src);
                    rhotracfab.resize(
                        rhotrac_box, PDE::ndim, amrex::The_Async_Arena());
                    rhotrac = rhotracfab.array();
                    rhotrac_nph_fab.resize(
                        rhotrac_box, PDE::ndim, amrex::The_Async_Arena());
                    rhotrac_nph = rhotrac_nph_fab.array();
                    if (mphase_vof) {
                        srctrac_over_rho_fab.resize(
                            rhotrac_box, PDE::ndim, amrex::The_Async_Arena());
                        srctrac_over_rho = srctrac_over_rho_fab.array();
                    }

                    amrex::ParallelFor(
                        rhotrac_box, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                            rhotrac(i, j, k, n) =
                                rho_arr(i, j, k) * tra_arr(i, j, k, n);
                            rhotrac_nph(i, j, k, n) =
                                rho_nph_arr(i, j, k) * tra_nph_arr(i, j, k, n);
                            if (mphase_vof && srctrac_box.contains(i, j, k)) {
                                // Need to remove the effect of the multiply_rho
                                // routine -- the setup is different for
                                // mass-consistent multiphase advection
                                srctrac_over_rho(i, j, k, n) =
                                    src_arr(i, j, k, n) / rho_nph_arr(i, j, k);
                            }
                        });
                }

                if ((godunov_scheme == godunov::scheme::PPM) ||
                    (godunov_scheme == godunov::scheme::PPM_NOLIM) ||
                    (godunov_scheme == godunov::scheme::PLM) ||
                    (godunov_scheme == godunov::scheme::WENOZ) ||
                    (godunov_scheme == godunov::scheme::WENO_JS) ||
                    (godunov_scheme == godunov::scheme::BDS)) {
                    amrex::FArrayBox tmpfab(
                        amrex::grow(bx, 1), 1, amrex::The_Async_Arena());
                    tmpfab.setVal<amrex::RunOn::Device>(0.0_rt);
                    const auto& divu = tmpfab.array();
                    const bool is_velocity = false;
                    const bool known_edge_state = false;
                    const bool godunov_use_ppm =
                        ((godunov_scheme != godunov::scheme::PLM) &&
                         (godunov_scheme != godunov::scheme::BDS));
                    int limiter_type;
                    if (godunov_scheme == godunov::scheme::PPM_NOLIM) {
                        limiter_type = PPM::NoLimiter;
                    } else if (godunov_scheme == godunov::scheme::WENOZ) {
                        limiter_type = PPM::WENOZ;
                    } else if (godunov_scheme == godunov::scheme::WENO_JS) {
                        limiter_type = PPM::WENO_JS;
                    } else {
                        limiter_type = PPM::default_limiter;
                    }

                    // if state is NPH, then n and n+1 are known, and only
                    // spatial extrapolation is performed
                    const amrex::Real dt_extrap =
                        (fstate == FieldState::NPH) ? 0.0_rt : dt;

                    HydroUtils::ComputeFluxesOnBoxFromState(
                        bx, PDE::ndim, mfi,
                        ((PDE::multiply_rho && !mphase_vof) ? rhotrac
                                                            : tra_arr),
                        ((PDE::multiply_rho && !mphase_vof) ? rhotrac_nph
                                                            : tra_nph_arr),
                        (*flux_x)(lev).array(mfi), (*flux_y)(lev).array(mfi),
                        (*flux_z)(lev).array(mfi), (*face_x)(lev).array(mfi),
                        (*face_y)(lev).array(mfi), (*face_z)(lev).array(mfi),
                        known_edge_state, u_mac(lev).const_array(mfi),
                        v_mac(lev).const_array(mfi),
                        w_mac(lev).const_array(mfi), divu,
                        ((PDE::multiply_rho && mphase_vof)
                             ? srctrac_over_rho
                             : src_term(lev).const_array(mfi)),
                        geom[lev], dt_extrap, dof_field.bcrec(),
                        dof_field.bcrec_device().data(), iconserv.data(),
                        godunov_use_ppm, godunov_use_forces_in_trans,
                        is_velocity, fluxes_are_area_weighted, advection_type,
                        limiter_type, m_allow_inflow_on_outflow);
                } else {
                    amrex::Abort("Invalid godunov scheme");
                }
            }
        }

        // Multiphase flux operations
        if (mphase_vof && PDE::multiply_rho) {
            // Loop levels
            multiphase::hybrid_fluxes(
                repo, PDE::ndim, iconserv, (*flux_x), (*flux_y), (*flux_z),
                dof_field, dof_nph, src_term, den, den_nph, u_mac, v_mac, w_mac,
                dof_field.bcrec(), dof_field.bcrec_device().data(), den.bcrec(),
                den.bcrec_device().data(), dt, mflux_scheme,
                m_allow_inflow_on_outflow, godunov_use_forces_in_trans, true);
        }

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
            repo.num_active_levels());
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            fluxes[lev][0] = &(*flux_x)(lev);
            fluxes[lev][1] = &(*flux_y)(lev);
            fluxes[lev][2] = &(*flux_z)(lev);
        }

        // In order to enforce conservation across coarse-fine boundaries we
        // must be sure to average down the fluxes before we use them
        for (int lev = repo.num_active_levels() - 1; lev > 0; --lev) {
            amrex::IntVect rr =
                geom[lev].Domain().size() / geom[lev - 1].Domain().size();
            amrex::average_down_faces(
                GetArrOfConstPtrs(fluxes[lev]), fluxes[lev - 1], rr,
                geom[lev - 1]);
        }

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();

                HydroUtils::ComputeDivergence(
                    bx, conv_term(lev).array(mfi), (*flux_x)(lev).array(mfi),
                    (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                    PDE::ndim, geom[lev], amrex::Real(-1.0_rt),
                    fluxes_are_area_weighted);
            }
        }
    }

    PDEFields& fields;
    Field& density;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    amrex::Gpu::DeviceVector<int> iconserv;

    godunov::scheme godunov_scheme = godunov::scheme::WENOZ;
    godunov::scheme mflux_scheme = godunov::scheme::UPWIND;
    std::string godunov_type{"weno_z"};
    std::string mflux_type{"upwind"};
    const bool fluxes_are_area_weighted{false};
    bool godunov_use_forces_in_trans{false};
    bool m_allow_inflow_on_outflow{false};
    std::string advection_type{"Godunov"};
};

} // namespace amr_wind::pde

#endif /* ADVOP_GODUNOV_H */
