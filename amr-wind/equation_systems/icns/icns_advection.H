#ifndef ICNS_ADVECTION_H
#define ICNS_ADVECTION_H

#include "amr-wind/equation_systems/AdvOp_Godunov.H"
#include "amr-wind/equation_systems/AdvOp_MOL.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/core/Physics.H"

#include "AMReX_MultiFabUtil.H"
#include "hydro_MacProjector.H"
#include "hydro_mol.H"
#include "hydro_utils.H"

#include "amr-wind/utilities/diagnostics.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

class MacProjOp
{
public:
    using FaceFabPtrVec =
        amrex::Vector<amrex::Array<const amrex::MultiFab*, ICNS::ndim>>;

    MacProjOp(
        FieldRepo& /*repo*/,
        PhysicsMgr& /*phy_mgr*/,
        bool /*has_overset*/,
        bool /*variable_density*/,
        bool /*mesh_mapping*/,
        bool /*is_anelastic*/);

    void set_inflow_velocity(amrex::Real time);

    void operator()(FieldState fstate, amrex::Real dt);

    static void mac_proj_to_uniform_space(
        const amr_wind::FieldRepo& /*repo*/,
        amr_wind::Field& /*u_mac*/,
        amr_wind::Field& /*v_mac*/,
        amr_wind::Field& /*w_mac*/,
        amrex::Array<amrex::MultiFab*, ICNS::ndim>& /*rho_face*/,
        amrex::Real /*ovst_fac*/,
        int /*lev*/);

    [[nodiscard]] amrex::Real rho0() const { return m_rho_0; }

private:
    void init_projector(const FaceFabPtrVec& /*beta*/);
    void init_projector(amrex::Real /*beta*/);

    void enforce_inout_solvability(
        const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
            a_umac);

    FieldRepo& m_repo;
    PhysicsMgr& m_phy_mgr;
    std::unique_ptr<Hydro::MacProjector> m_mac_proj;
#ifdef AMR_WIND_USE_FFT
    std::unique_ptr<Hydro::FFTMacProjector> m_fft_mac_proj;
    bool m_use_fft{true}; // use fft if possible
#endif
    MLMGOptions m_options;
    bool m_has_overset{false};
    bool m_need_init{true};
    bool m_variable_density{false};
    bool m_mesh_mapping{false};
    bool m_is_anelastic{false};
    amrex::Real m_rho_0{1.0_rt};
};

/** Godunov scheme for ICNS
 *  \ingroup icns
 */
template <>
struct AdvectionOp<ICNS, fvm::Godunov>
{
    AdvectionOp(
        CFDSim& sim,
        PDEFields& fields_in,
        bool has_overset,
        bool variable_density,
        bool mesh_mapping,
        bool is_anelastic)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
        , m_macproj_op(
              fields.repo,
              sim.physics_manager(),
              has_overset,
              variable_density,
              mesh_mapping,
              is_anelastic)
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
            // use Godunov for premac, use BDS for postmac. Eventually
            // there will be a premac BDS
            postmac_advection_type = "BDS";
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

        // Formulation of discrete ICNS equation
        // 1 = conservative (default), 0 = nonconservative
        pp.query("icns_conserv", m_cons);
        iconserv.resize(ICNS::ndim, m_cons);

        // Get copy of verbose
        pp.query("verbose", m_verbose);

        amrex::ParmParse pp_eq("ICNS");
        pp_eq.query(
            "allow_inflow_at_pressure_outflow", m_allow_inflow_on_outflow);
    }

    void preadvect(
        const FieldState fstate, const amrex::Real dt, const amrex::Real time)
    {

        const auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        const auto& dof_field = fields.field.state(fstate);
        auto bcrec_device = dof_field.bcrec_device();

        //
        // Predict
        //
        if ((godunov_scheme == godunov::scheme::PPM) ||
            (godunov_scheme == godunov::scheme::PPM_NOLIM) ||
            (godunov_scheme == godunov::scheme::PLM) ||
            (godunov_scheme == godunov::scheme::WENOZ) ||
            (godunov_scheme == godunov::scheme::WENO_JS) ||
            (godunov_scheme == godunov::scheme::BDS)) {
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
            for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
                HydroUtils::ExtrapVelToFaces(
                    dof_field(lev), src_term(lev), u_mac(lev), v_mac(lev),
                    w_mac(lev), dof_field.bcrec(), bcrec_device.data(),
                    repo.mesh().Geom(lev), dt_extrap, godunov_use_ppm,
                    godunov_use_forces_in_trans, premac_advection_type,
                    limiter_type, m_allow_inflow_on_outflow);
            }
        } else {
            amrex::Abort("Invalid godunov scheme");
        }

        if (m_verbose > 2) {
            diagnostics::PrintMaxMACVelLocations(repo, "before MAC projection");
        }

        // Populate boundaries (valid cells) using velocity BCs
        m_macproj_op.set_inflow_velocity(time);

        // MAC projection
        m_macproj_op(fstate, dt);

        // Fill mac velocities (ghost cells) using velocity BCs
        if (fvm::Godunov::nghost_state > 0) {
            amrex::Array<Field*, AMREX_SPACEDIM> mac_vel = {
                AMREX_D_DECL(&u_mac, &v_mac, &w_mac)};
            dof_field.fillpatch_sibling_fields(time, u_mac.num_grow(), mac_vel);
        }

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            u_mac(lev).FillBoundary(geom[lev].periodicity());
            v_mac(lev).FillBoundary(geom[lev].periodicity());
            w_mac(lev).FillBoundary(geom[lev].periodicity());
        }

        if (m_verbose > 2) {
            diagnostics::PrintMaxMACVelLocations(repo, "after MAC projection");
        }
    }

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        const auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        const auto& src_term = fields.src_term;
        // cppcheck-suppress constVariableReference
        auto& conv_term = fields.conv_term;
        const auto& dof_field = fields.field.state(fstate);
        const auto& dof_nph = fields.field.state(amr_wind::FieldState::NPH);

        auto flux_x =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto flux_y =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto flux_z =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::ZFACE);
        auto face_x =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::XFACE);
        auto face_y =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::YFACE);
        auto face_z =
            repo.create_scratch_field(ICNS::ndim, 0, amr_wind::FieldLoc::ZFACE);

        const auto& rho_o =
            repo.get_field("density").state(amr_wind::FieldState::Old);
        const auto& rho_nph =
            repo.get_field("density").state(amr_wind::FieldState::NPH);

        const bool mphase_vof = repo.field_exists("vof");

        //
        // Advect momentum eqns
        //
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {

            // form multifab for transport variable and source term
            amrex::MultiFab q(
                dof_field(lev).boxArray(), dof_field(lev).DistributionMap(),
                ICNS::ndim, fvm::Godunov::nghost_state);
            amrex::MultiFab::Copy(
                q, dof_field(lev), 0, 0, ICNS::ndim,
                fvm::Godunov::nghost_state);
            amrex::MultiFab fq(
                src_term(lev).boxArray(), src_term(lev).DistributionMap(),
                ICNS::ndim, fvm::Godunov::nghost_src);
            amrex::MultiFab::Copy(
                fq, src_term(lev), 0, 0, ICNS::ndim, fvm::Godunov::nghost_src);
            // form multifab for time-correct boundary condition of variable
            amrex::MultiFab q_nph(
                dof_field(lev).boxArray(), dof_field(lev).DistributionMap(),
                ICNS::ndim, fvm::Godunov::nghost_state);
            amrex::MultiFab::Copy(
                q_nph, dof_nph(lev), 0, 0, ICNS::ndim,
                fvm::Godunov::nghost_state);

            // Calculate fluxes using momentum directly
            if (!mphase_vof) {

                for (int idim = 0; idim < dof_field.num_comp(); ++idim) {
                    amrex::MultiFab::Multiply(
                        q, rho_o(lev), 0, idim, 1, fvm::Godunov::nghost_state);
                    // Source terms at old state during advection calculation
                    amrex::MultiFab::Multiply(
                        fq, rho_o(lev), 0, idim, 1, fvm::Godunov::nghost_src);

                    amrex::MultiFab::Multiply(
                        q_nph, rho_nph(lev), 0, idim, 1,
                        fvm::Godunov::nghost_state);
                }
            }

            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
            }
            if ((godunov_scheme == godunov::scheme::PPM) ||
                (godunov_scheme == godunov::scheme::PPM_NOLIM) ||
                (godunov_scheme == godunov::scheme::PLM) ||
                (godunov_scheme == godunov::scheme::WENOZ) ||
                (godunov_scheme == godunov::scheme::WENO_JS) ||
                (godunov_scheme == godunov::scheme::BDS)) {
                const bool is_velocity = true;
                const bool known_edge_state = false;
                const bool godunov_use_ppm =
                    ((godunov_scheme == godunov::scheme::PPM) ||
                     (godunov_scheme == godunov::scheme::PPM_NOLIM) ||
                     (godunov_scheme == godunov::scheme::WENOZ) ||
                     (godunov_scheme == godunov::scheme::WENO_JS));
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
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                     ++mfi) {
                    const auto& bx = mfi.tilebox();
                    amrex::FArrayBox tmpfab(
                        amrex::grow(bx, 1), 1, amrex::The_Async_Arena());
                    tmpfab.setVal<amrex::RunOn::Device>(0.0_rt);
                    const auto& divu = tmpfab.array();
                    HydroUtils::ComputeFluxesOnBoxFromState(
                        bx, ICNS::ndim, mfi, q.const_array(mfi),
                        q_nph.const_array(mfi), (*flux_x)(lev).array(mfi),
                        (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                        (*face_x)(lev).array(mfi), (*face_y)(lev).array(mfi),
                        (*face_z)(lev).array(mfi), known_edge_state,
                        u_mac(lev).const_array(mfi),
                        v_mac(lev).const_array(mfi),
                        w_mac(lev).const_array(mfi), divu, fq.const_array(mfi),
                        geom[lev], dt_extrap, dof_field.bcrec(),
                        dof_field.bcrec_device().data(), iconserv.data(),
                        godunov_use_ppm, godunov_use_forces_in_trans,
                        is_velocity, fluxes_are_area_weighted,
                        postmac_advection_type, limiter_type,
                        m_allow_inflow_on_outflow);
                }
            } else {
                amrex::Abort("Invalid godunov scheme");
            }
        }

        // Multiphase flux operations
        if (mphase_vof) {
            // Loop levels
            multiphase::hybrid_fluxes(
                repo, ICNS::ndim, iconserv, (*flux_x), (*flux_y), (*flux_z),
                dof_field, dof_nph, src_term, rho_o, rho_nph, u_mac, v_mac,
                w_mac, dof_field.bcrec(), dof_field.bcrec_device().data(),
                rho_o.bcrec(), rho_o.bcrec_device().data(), dt, mflux_scheme,
                m_allow_inflow_on_outflow, godunov_use_forces_in_trans);
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
                    ICNS::ndim, geom[lev], static_cast<amrex::Real>(-1.0_rt),
                    fluxes_are_area_weighted);

                if (m_cons == 0) {
                    amrex::FArrayBox div_umac(bx, 1, amrex::The_Async_Arena());
                    auto const& divum_arr = div_umac.array();
                    HydroUtils::ComputeDivergence(
                        bx, divum_arr, u_mac(lev).const_array(mfi),
                        v_mac(lev).const_array(mfi),
                        w_mac(lev).const_array(mfi), 1, geom[lev],
                        static_cast<amrex::Real>(1.0_rt), false);
                    HydroUtils::ComputeConvectiveTerm(
                        bx, ICNS::ndim, mfi, dof_field(lev).const_array(mfi),
                        (*face_x)(lev).const_array(mfi),
                        (*face_y)(lev).const_array(mfi),
                        (*face_z)(lev).const_array(mfi), divum_arr,
                        conv_term(lev).array(mfi), iconserv.data(),
                        postmac_advection_type);
                }
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;

    MacProjOp m_macproj_op;
    amrex::Gpu::DeviceVector<int> iconserv;

    godunov::scheme godunov_scheme = godunov::scheme::WENOZ;
    godunov::scheme mflux_scheme = godunov::scheme::UPWIND;
    std::string godunov_type{"weno_z"};
    std::string mflux_type{"upwind"};
    const bool fluxes_are_area_weighted{false};
    bool godunov_use_forces_in_trans{false};
    int m_cons{1};
    int m_verbose{0};
    bool m_allow_inflow_on_outflow{false};
    std::string premac_advection_type{"Godunov"};
    std::string postmac_advection_type{"Godunov"};
};

/** MOL scheme for ICNS
 *  \ingroup icns
 */
template <>
struct AdvectionOp<ICNS, fvm::MOL>
{
    AdvectionOp(
        CFDSim& sim,
        PDEFields& fields_in,
        bool has_overset,
        bool variable_density,
        bool mesh_mapping,
        bool is_anelastic)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
        , m_mesh_mapping(mesh_mapping)
        , m_macproj_op(
              fields.repo,
              sim.physics_manager(),
              has_overset,
              variable_density,
              m_mesh_mapping,
              is_anelastic)
    {}

    void preadvect(
        const FieldState fstate,
        const amrex::Real dt,
        const amrex::Real /*time*/)
    {

        const auto& repo = fields.repo;
        auto& dof_field = fields.field.state(fstate);

        // computation of velocity on faces requires
        // dof field to be in stretched mesh space
        if (dof_field.in_uniform_space() && m_mesh_mapping) {
            dof_field.to_stretched_space();
        }

        //
        // Predict velocities
        //

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            MOL::ExtrapVelToFaces(
                dof_field(lev), u_mac(lev), v_mac(lev), w_mac(lev),
                repo.mesh().Geom(lev), dof_field.bcrec(),
                dof_field.bcrec_device().data());
        }

        m_macproj_op(fstate, dt);
    }

    void operator()(const FieldState fstate, const amrex::Real /*unused*/)
    {

        const auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();
        // cppcheck-suppress constVariableReference
        auto& conv_term = fields.conv_term.state(fstate);
        const auto& dof_field = fields.field.state(fstate);
        const auto& rho = repo.get_field("density").state(fstate);

        //
        // Advect velocity
        //

        int nmaxcomp = AMREX_SPACEDIM;
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {

            amrex::MFItInfo mfi_info;
            // if (amrex::Gpu::notInLaunchRegion())
            // mfi_info.EnableTiling(amrex::IntVect(1024,16,16)).SetDynamic(true);
            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                amrex::Box const& bx = mfi.tilebox();
                amrex::Box gbx = grow(bx, fvm::MOL::nghost_state);

                // Set up momentum array
                amrex::FArrayBox qfab(
                    gbx, ICNS::ndim, amrex::The_Async_Arena());
                const auto& q = qfab.array();
                // Calculate momentum
                auto rho_arr = rho(lev).const_array(mfi);
                auto vel_arr = dof_field(lev).const_array(mfi);
                amrex::ParallelFor(
                    gbx, ICNS::ndim,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                        q(i, j, k, n) = rho_arr(i, j, k) * vel_arr(i, j, k, n);
                    });
                // Doing this explicitly, instead of through a Multiply command,
                // helps avoid floating-point errors with intel compilers and
                // mimics the implementation in equation_systems/AdvOp_MOL.H

                amrex::Box tmpbox = amrex::surroundingNodes(bx);
                const int tmpcomp = nmaxcomp * AMREX_SPACEDIM;

                amrex::FArrayBox tmpfab(
                    tmpbox, tmpcomp, amrex::The_Async_Arena());

                amrex::Array4<amrex::Real> fx = tmpfab.array(0);
                amrex::Array4<amrex::Real> fy = tmpfab.array(nmaxcomp);
                amrex::Array4<amrex::Real> fz = tmpfab.array(nmaxcomp * 2);

                mol::compute_convective_fluxes(
                    lev, bx, AMREX_SPACEDIM, fx, fy, fz, q,
                    u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                    w_mac(lev).const_array(mfi), dof_field.bcrec().data(),
                    dof_field.bcrec_device().data(), geom);

                mol::compute_convective_rate(
                    bx, AMREX_SPACEDIM, conv_term(lev).array(mfi), fx, fy, fz,
                    geom[lev].InvCellSizeArray());
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;

    bool m_mesh_mapping;

    MacProjOp m_macproj_op;
};

} // namespace amr_wind::pde

#endif /* ICNS_ADVECTION_H */
