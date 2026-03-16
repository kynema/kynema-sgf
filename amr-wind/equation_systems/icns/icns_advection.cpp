#include <memory>
#include "AMReX.H"
#include "amr-wind/equation_systems/icns/icns_advection.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/overset/overset_ops_routines.H"
#include "AMReX_MultiFabUtil.H"
#include "hydro_MacProjector.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

namespace {

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> get_projection_bc(
    amrex::Orientation::Side side,
    amrex::GpuArray<BC, AMREX_SPACEDIM * 2> bctype,
    amrex::Vector<amrex::Geometry> geom)
{

    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = amrex::LinOpBCType::Periodic;
        } else {
            auto bc = bctype[amrex::Orientation(dir, side)];
            if (bc == BC::pressure_outflow) {
                r[dir] = amrex::LinOpBCType::Dirichlet;
            } else {
                r[dir] = amrex::LinOpBCType::Neumann;
            }
        }
    }
    return r;
}

} // namespace

MacProjOp::MacProjOp(
    FieldRepo& repo,
    PhysicsMgr& phy_mgr,
    bool has_overset,
    bool variable_density,
    bool mesh_mapping,
    bool is_anelastic)
    : m_repo(repo)
    , m_phy_mgr(phy_mgr)
    , m_options("mac_proj")
    , m_has_overset(has_overset)
    , m_variable_density(variable_density)
    , m_mesh_mapping(mesh_mapping)
    , m_is_anelastic(is_anelastic)
{
    amrex::ParmParse pp("incflo");
    pp.query("density", m_rho_0);
#ifdef AMR_WIND_USE_FFT
    {
        amrex::ParmParse pp("mac_proj");
        pp.query("use_fft", m_use_fft);
    }
#endif
}

void MacProjOp::enforce_inout_solvability(
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>& a_umac)
{
    auto& velocity = m_repo.get_field("velocity");
    amrex::BCRec const* bc_type = velocity.bcrec().data();
    const amrex::Vector<amrex::Geometry>& geom = m_repo.mesh().Geom();

    HydroUtils::enforceInOutSolvability(a_umac, bc_type, geom);
}

void MacProjOp::init_projector(const MacProjOp::FaceFabPtrVec& beta)
{
    // Prepare masking for projection
    if (m_has_overset) {
        amr_wind::overset_ops::prepare_mask_cell_for_mac(m_repo);
    }

    // Prepare projector
    m_mac_proj = std::make_unique<Hydro::MacProjector>(
        m_repo.mesh().Geom(0, m_repo.num_active_levels() - 1));
    m_mac_proj->initProjector(
        m_options.lpinfo(), beta,
        m_has_overset ? m_repo.get_int_field("mask_cell").vec_const_ptrs()
                      : amrex::Vector<const amrex::iMultiFab*>());

    m_options(*m_mac_proj);

    auto& pressure = m_repo.get_field("p");
    const auto& bctype = pressure.bc_type();

    m_mac_proj->setDomainBC(
        get_projection_bc(
            amrex::Orientation::low, bctype, m_repo.mesh().Geom()),
        get_projection_bc(
            amrex::Orientation::high, bctype, m_repo.mesh().Geom()));

    m_need_init = false;
}

void MacProjOp::init_projector(const amrex::Real beta)
{
    // Prepare masking for projection
    if (m_has_overset) {
        amr_wind::overset_ops::prepare_mask_cell_for_mac(m_repo);
    }

    // Prepare projector
    m_need_init = false;

    auto& pressure = m_repo.get_field("p");
    const auto& bctype = pressure.bc_type();

    auto const& lobc = get_projection_bc(
        amrex::Orientation::low, bctype, m_repo.mesh().Geom());
    auto const& hibc = get_projection_bc(
        amrex::Orientation::high, bctype, m_repo.mesh().Geom());

#ifdef AMR_WIND_USE_FFT
    if (m_use_fft) {
        if (m_repo.num_active_levels() == 1 && m_has_overset == false) {
            m_fft_mac_proj = std::make_unique<Hydro::FFTMacProjector>(
                m_repo.mesh().Geom(0), lobc, hibc);
            return;
        } else {
            amrex::ParmParse pp("mac_proj");
            if (pp.contains("use_fft")) {
                amrex::Print() << "WARNING: FFT MAC projection disabled due to "
                                  "multiple levels/overset\n";
            }
        }
    }
#endif

    m_mac_proj = std::make_unique<Hydro::MacProjector>(
        m_repo.mesh().Geom(0, m_repo.num_active_levels() - 1));
    m_mac_proj->initProjector(
        m_repo.mesh().boxArray(0, m_repo.num_active_levels() - 1),
        m_repo.mesh().DistributionMap(0, m_repo.num_active_levels() - 1),
        m_options.lpinfo(), beta,
        m_has_overset ? m_repo.get_int_field("mask_cell").vec_const_ptrs()
                      : amrex::Vector<const amrex::iMultiFab*>());

    m_options(*m_mac_proj);

    m_mac_proj->setDomainBC(lobc, hibc);
}

void MacProjOp::set_inflow_velocity(amrex::Real time)
{
    // Currently, input boundary planes account for inflow differently
    // Also, MPL needs to be refactored to do this properly, defer for now
    if (m_phy_mgr.contains("ABL")) {
        if (m_phy_mgr.get<amr_wind::ABL>().bndry_plane().mode() ==
            io_mode::input) {
            return;
        }
        if (m_phy_mgr.get<amr_wind::ABL>().abl_mpl().is_active()) {
            return;
        }
    }

    auto& velocity = m_repo.get_field("velocity");
    auto& u_mac = m_repo.get_field("u_mac");
    auto& v_mac = m_repo.get_field("v_mac");
    auto& w_mac = m_repo.get_field("w_mac");

    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> mac_vec = {
            AMREX_D_DECL(&u_mac(lev), &v_mac(lev), &w_mac(lev))};
        velocity.set_inflow_sibling_fields(lev, time, mac_vec);
        if (m_phy_mgr.contains("OceanWaves")) {
            auto& ow = m_phy_mgr.get<amr_wind::ocean_waves::OceanWaves>();
            ow.ow_bndry().set_inflow_sibling_velocity(
                lev, time, velocity, mac_vec);
        }
    }
}

//
// Computes the following decomposition:
//
//    u + c*grad(phi)/ro = u*  with  div(ep*u) = 0
//
// Inputs:
//
//   u_mac,v_mac,w_mac = the MAC velocity field to be projected
//   density           = the cell-centered density
//
// Outputs:
//
//  u_mac,v_mac,w_mac = the PROJECTED MAC velocity field
//
// Notes:
//
//  phi, the projection auxiliary function, is computed by solving
//
//       div(ep*grad(phi)/rho) = div(ep * u*)
//

void MacProjOp::operator()(const FieldState fstate, const amrex::Real dt)
{
    BL_PROFILE("amr-wind::ICNS::advection_mac_project");
    const auto& geom = m_repo.mesh().Geom();
    auto& u_mac = m_repo.get_field("u_mac");
    auto& v_mac = m_repo.get_field("v_mac");
    auto& w_mac = m_repo.get_field("w_mac");
    const auto& density = m_repo.get_field("density", fstate);

    amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> mac_vec(
        m_repo.num_active_levels());

    amrex::Real factor = m_has_overset ? 0.5_rt * dt : 1.0_rt;

    std::unique_ptr<ScratchField> ref_rho_xf, ref_rho_yf, ref_rho_zf;
    if (m_is_anelastic) {
        ref_rho_xf =
            m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
        ref_rho_yf =
            m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
        ref_rho_zf =
            m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);
    }
    amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> ref_rho_face(
        m_repo.num_active_levels());

    // TODO: remove the or in the if statement for m_has_overset
    // For now assume variable viscosity for overset
    // this can be removed once the nsolve overset
    // masking is implemented in cell based AMReX poisson solvers
    if (m_variable_density || m_has_overset || m_mesh_mapping ||
        m_is_anelastic) {
        amrex::Vector<amrex::Array<amrex::MultiFab const*, ICNS::ndim>>
            rho_face_const;
        rho_face_const.reserve(m_repo.num_active_levels());

        // This will hold density on faces
        std::unique_ptr<ScratchField> rho_xf, rho_yf, rho_zf;
        rho_xf = m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
        rho_yf = m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
        rho_zf = m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);
        amrex::Vector<amrex::Array<amrex::MultiFab*, ICNS::ndim>> rho_face(
            m_repo.num_active_levels());

        const auto* ref_density =
            m_is_anelastic ? &(m_repo.get_field("reference_density")) : nullptr;

        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            rho_face[lev][0] = &(*rho_xf)(lev);
            rho_face[lev][1] = &(*rho_yf)(lev);
            rho_face[lev][2] = &(*rho_zf)(lev);

            amrex::average_cellcenter_to_face(
                rho_face[lev], density(lev), geom[lev]);

            if (m_is_anelastic) {
                ref_rho_face[lev][0] = &(*ref_rho_xf)(lev);
                ref_rho_face[lev][1] = &(*ref_rho_yf)(lev);
                ref_rho_face[lev][2] = &(*ref_rho_zf)(lev);
                amrex::average_cellcenter_to_face(
                    ref_rho_face[lev], (*ref_density)(lev), geom[lev]);
                for (int idim = 0; idim < ICNS::ndim; ++idim) {
                    rho_face[lev][idim]->divide(
                        *(ref_rho_face[lev][idim]), 0, density.num_comp(),
                        rho_face[lev][idim]->nGrow());
                }
            }

            if (m_mesh_mapping) {
                mac_proj_to_uniform_space(
                    m_repo, u_mac, v_mac, w_mac, rho_face[lev], factor, lev);
            } else {
                for (int idim = 0; idim < ICNS::ndim; ++idim) {
                    rho_face[lev][idim]->invert(factor, 0);
                }
            }

            rho_face_const.push_back(GetArrOfConstPtrs(rho_face[lev]));
        }

        if (m_need_init) {
            init_projector(rho_face_const);
        } else {
            m_mac_proj->updateBeta(rho_face_const);
        }
    } else {
        if (m_need_init) {
            init_projector(factor / m_rho_0);
        } else {
#ifdef AMR_WIND_USE_FFT
            if (!m_fft_mac_proj)
#endif
            {
                m_mac_proj->updateBeta(factor / m_rho_0);
            }
        }
    }

    for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
        mac_vec[lev][0] = &u_mac(lev);
        mac_vec[lev][1] = &v_mac(lev);
        mac_vec[lev][2] = &w_mac(lev);
        if (m_is_anelastic) {
            for (int idim = 0; idim < ICNS::ndim; ++idim) {
                amrex::Multiply(
                    *(mac_vec[lev][idim]), *(ref_rho_face[lev][idim]), 0, 0,
                    density.num_comp(), 0);
            }
        }
    }

    const bool has_inout_bndry =
        (m_repo.get_field("velocity")).has_inout_bndry();
    if (has_inout_bndry) {
        enforce_inout_solvability(mac_vec);
    }

#ifdef AMR_WIND_USE_FFT
    if (m_fft_mac_proj) {
        // This is set on mac_vec[0] since FFT based projection is restricted to
        // a single level
        m_fft_mac_proj->setUMAC(mac_vec[0]);
    } else
#endif
    {
        m_mac_proj->setUMAC(mac_vec);
    }

    if (m_has_overset) {
        // In masked regions, the pressure should not change from what was used
        // when preparing the MAC velocity field; therefore, phi is set to 0
        auto phif = m_repo.create_scratch_field(1, 1, amr_wind::FieldLoc::CELL);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*phif)(lev).setVal(0.);
        }

        m_mac_proj->project(
            phif->vec_ptrs(), m_options.rel_tol, m_options.abs_tol);

    } else {
#ifdef AMR_WIND_USE_FFT
        if (m_fft_mac_proj) {
            m_fft_mac_proj->project();
        } else
#endif
        {
            m_mac_proj->project(m_options.rel_tol, m_options.abs_tol);
        }
    }

    if (m_is_anelastic) {
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            for (int idim = 0; idim < ICNS::ndim; ++idim) {
                amrex::Divide(
                    *(mac_vec[lev][idim]), *(ref_rho_face[lev][idim]), 0, 0,
                    density.num_comp(), 0);
            }
        }
    }

    // Prepare masking for remainder of algorithm
    if (m_has_overset) {
        amr_wind::overset_ops::revert_mask_cell_after_mac(m_repo);
    }

    if (m_mac_proj) {
        io::print_mlmg_info("MAC_projection", m_mac_proj->getMLMG());
    }
}

void MacProjOp::mac_proj_to_uniform_space(
    const amr_wind::FieldRepo& repo,
    amr_wind::Field& u_mac,
    amr_wind::Field& v_mac,
    amr_wind::Field& w_mac,
    amrex::Array<amrex::MultiFab*, ICNS::ndim>& rho_face,
    amrex::Real ovst_fac,
    int lev)
{
    const auto& mesh_fac_xf =
        repo.get_mesh_mapping_field(amr_wind::FieldLoc::XFACE);
    const auto& mesh_fac_yf =
        repo.get_mesh_mapping_field(amr_wind::FieldLoc::YFACE);
    const auto& mesh_fac_zf =
        repo.get_mesh_mapping_field(amr_wind::FieldLoc::ZFACE);
    const auto& mesh_detJ_xf =
        repo.get_mesh_mapping_det_j(amr_wind::FieldLoc::XFACE);
    const auto& mesh_detJ_yf =
        repo.get_mesh_mapping_det_j(amr_wind::FieldLoc::YFACE);
    const auto& mesh_detJ_zf =
        repo.get_mesh_mapping_det_j(amr_wind::FieldLoc::ZFACE);

    // scale U^mac to accommodate for mesh mapping -> U^bar = J/fac *
    // U^mac beta accounted for mesh mapping = J/fac^2 * 1/rho construct
    // rho and mesh map u_mac on x-face
    {
        const auto& u_arrs = u_mac(lev).arrays();
        const auto& rho_arrs = rho_face[0]->arrays();
        const auto& fac_arrs = mesh_fac_xf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_xf(lev).const_arrays();

        amrex::ParallelFor(
            *(rho_face[0]), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                u_arrs[nbx](i, j, k) *=
                    detJ_arrs[nbx](i, j, k) / fac_arrs[nbx](i, j, k, 0);
                rho_arrs[nbx](i, j, k) =
                    ovst_fac * detJ_arrs[nbx](i, j, k) /
                    std::pow(fac_arrs[nbx](i, j, k, 0), 2.0_rt) /
                    rho_arrs[nbx](i, j, k);
            });
    }

    // construct rho on y-face
    {
        const auto& v_arrs = v_mac(lev).arrays();
        const auto& rho_arrs = rho_face[1]->arrays();
        const auto& fac_arrs = mesh_fac_yf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_yf(lev).const_arrays();

        amrex::ParallelFor(
            *(rho_face[1]), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                v_arrs[nbx](i, j, k) *=
                    detJ_arrs[nbx](i, j, k) / fac_arrs[nbx](i, j, k, 1);
                rho_arrs[nbx](i, j, k) =
                    ovst_fac * detJ_arrs[nbx](i, j, k) /
                    std::pow(fac_arrs[nbx](i, j, k, 1), 2.0_rt) /
                    rho_arrs[nbx](i, j, k);
            });
    }

    // construct rho on z-face
    {
        const auto& w_arrs = w_mac(lev).arrays();
        const auto& rho_arrs = rho_face[2]->arrays();
        const auto& fac_arrs = mesh_fac_zf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_zf(lev).const_arrays();

        amrex::ParallelFor(
            *(rho_face[2]), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                w_arrs[nbx](i, j, k) *=
                    detJ_arrs[nbx](i, j, k) / fac_arrs[nbx](i, j, k, 2);
                rho_arrs[nbx](i, j, k) =
                    ovst_fac * detJ_arrs[nbx](i, j, k) /
                    std::pow(fac_arrs[nbx](i, j, k, 2), 2.0_rt) /
                    rho_arrs[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind::pde
