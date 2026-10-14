#include "src/incflo.H"
#include "src/diffusion/diffusion.H"
#include "AMReX_REAL.H"
#include "src/utilities/math_ops.H"

using namespace amrex::literals;

namespace diffusion {

amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
get_diffuse_tensor_bc(
    kynema_sgf::Field& velocity, amrex::Orientation::Side side)
{
    const auto& geom = velocity.repo().mesh().Geom(0);
    amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>> r(3);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom.isPeriodic(dir)) {
            r[0][dir] = amrex::LinOpBCType::Periodic;
            r[1][dir] = amrex::LinOpBCType::Periodic;
            r[2][dir] = amrex::LinOpBCType::Periodic;
        } else {
            auto bc = velocity.bc_type()[amrex::Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_outflow:
            case BC::zero_gradient: {
                // All three components are Neumann
                r[0][dir] = amrex::LinOpBCType::Neumann;
                r[1][dir] = amrex::LinOpBCType::Neumann;
                r[2][dir] = amrex::LinOpBCType::Neumann;
                break;
            }
            case BC::wave_generation:
            case BC::mass_inflow:
            case BC::mass_inflow_outflow:
            case BC::no_slip_wall: {
                // All three components are Dirichlet
                r[0][dir] = amrex::LinOpBCType::Dirichlet;
                r[1][dir] = amrex::LinOpBCType::Dirichlet;
                r[2][dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            case BC::symmetric_wall:
            case BC::slip_wall: {
                // Tangential components are Neumann
                // Normal     component  is  Dirichlet
                r[0][dir] = amrex::LinOpBCType::Neumann;
                r[1][dir] = amrex::LinOpBCType::Neumann;
                r[2][dir] = amrex::LinOpBCType::Neumann;

                r[dir][dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            case BC::wall_model: {
                // Tangential components are inhomogeneous Neumann
                // Normal     component  is  Dirichlet
                r[0][dir] = amrex::LinOpBCType::inhomogNeumann;
                r[1][dir] = amrex::LinOpBCType::inhomogNeumann;
                r[2][dir] = amrex::LinOpBCType::inhomogNeumann;

                r[dir][dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return r;
}

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>
get_diffuse_scalar_bc(kynema_sgf::Field& scalar, amrex::Orientation::Side side)
{
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (scalar.repo().mesh().Geom(0).isPeriodic(dir)) {
            r[dir] = amrex::LinOpBCType::Periodic;
        } else {
            auto bc = scalar.bc_type()[amrex::Orientation(dir, side)];
            switch (bc) {
            case BC::pressure_outflow:
            case BC::zero_gradient:
            case BC::symmetric_wall:
            case BC::slip_wall: {
                r[dir] = amrex::LinOpBCType::Neumann;
                break;
            }
            case BC::wall_model:
            case BC::fixed_gradient: {
                r[dir] = amrex::LinOpBCType::inhomogNeumann;
                break;
            }
            case BC::wave_generation:
            case BC::mass_inflow:
            case BC::mass_inflow_outflow:
            case BC::no_slip_wall: {
                r[dir] = amrex::LinOpBCType::Dirichlet;
                break;
            }
            default:
                amrex::Abort("get_diffuse_tensor_bc: undefined BC type");
            };
        }
    }
    return r;
}

amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> average_velocity_eta_to_faces(
    const amrex::Geometry& geom, amrex::MultiFab const& cc_eta)
{
    BL_PROFILE("kynema-sgf::diffusion::average_velocity_eta_to_faces");
    const auto& ba = cc_eta.boxArray();
    const auto& dm = cc_eta.DistributionMap();
    const auto& fact = cc_eta.Factory();
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> r{
        {amrex::MultiFab(
             amrex::convert(ba, amrex::IntVect::TheDimensionVector(0)), dm, 1,
             0, amrex::MFInfo(), fact),
         amrex::MultiFab(
             amrex::convert(ba, amrex::IntVect::TheDimensionVector(1)), dm, 1,
             0, amrex::MFInfo(), fact),
         amrex::MultiFab(
             amrex::convert(ba, amrex::IntVect::TheDimensionVector(2)), dm, 1,
             0, amrex::MFInfo(), fact)}};
    amrex::average_cellcenter_to_face(GetArrOfPtrs(r), cc_eta, geom);
    fixup_eta_on_domain_faces(geom, r, cc_eta);
    return r;
}

void fixup_eta_on_domain_faces(
    const amrex::Geometry& geom,
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fc,
    amrex::MultiFab const& cc)
{
    BL_PROFILE("kynema-sgf::diffusion::fixup_eta_on_domain_faces");

    const amrex::Box& domain = geom.Domain();
    amrex::MFItInfo mfi_info{};
    if (amrex::Gpu::notInLaunchRegion()) {
        mfi_info.SetDynamic(true);
    }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (false)
#endif
    for (amrex::MFIter mfi(cc, mfi_info); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.validbox();
        const auto& cca = cc.const_array(mfi);

        int idim = 0;
        if (!geom.isPeriodic(idim)) {
            const auto& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        fca(i, j, k) = cca(i, j, k);
                    });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        fca(i, j, k) = cca(i - 1, j, k);
                    });
            }
        }

        idim = 1;
        if (!geom.isPeriodic(idim)) {
            const auto& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        fca(i, j, k) = cca(i, j, k);
                    });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        fca(i, j, k) = cca(i, j - 1, k);
                    });
            }
        }

        idim = 2;
        if (!geom.isPeriodic(idim)) {
            const auto& fca = fc[idim].array(mfi);
            if (bx.smallEnd(idim) == domain.smallEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryLo(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        fca(i, j, k) = cca(i, j, k);
                    });
            }
            if (bx.bigEnd(idim) == domain.bigEnd(idim)) {
                amrex::ParallelFor(
                    amrex::bdryHi(bx, idim),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        fca(i, j, k) = cca(i, j, k - 1);
                    });
            }
        }
    }
}

void apply_immersed_interface(
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& fc,
    const kynema_sgf::FieldRepo& repo,
    const int lev)
{
    if (!repo.field_exists("terrain_diffusion_xf")) {
        return;
    }
    BL_PROFILE("kynema-sgf::diffusion::apply_immersed_interface");
    const amrex::Array<std::string, AMREX_SPACEDIM> names{
        {"terrain_diffusion_xf", "terrain_diffusion_yf",
         "terrain_diffusion_zf"}};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::MultiFab::Multiply(
            fc[dir], repo.get_field(names[dir])(lev), 0, 0, 1, 0);
    }
}

std::unique_ptr<kynema_sgf::ScratchField> immersed_effective_density(
    const kynema_sgf::FieldRepo& repo,
    const kynema_sgf::Field& density,
    const amrex::Real dt)
{
    if (!repo.field_exists("terrain_drag_rate")) {
        return nullptr;
    }
    BL_PROFILE("kynema-sgf::diffusion::immersed_effective_density");
    const auto& drag_rate = repo.get_field("terrain_drag_rate");
    auto rho_eff = repo.create_scratch_field(1, density.num_grow()[0]);
    // The drag-rate field carries fewer ghost cells than the density; fill
    // the ghosts with the plain density and scale only where the rate exists
    const amrex::IntVect ngrow_rate =
        amrex::min(density.num_grow(), drag_rate.num_grow());
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Copy(
            (*rho_eff)(lev), density(lev), 0, 0, 1, density.num_grow());
        const auto& rate_arrs = drag_rate(lev).const_arrays();
        const auto& reff_arrs = (*rho_eff)(lev).arrays();
        amrex::ParallelFor(
            (*rho_eff)(lev), ngrow_rate,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                reff_arrs[nbx](i, j, k) *=
                    (1.0_rt + (dt * rate_arrs[nbx](i, j, k)));
            });
    }
    amrex::Gpu::streamSynchronize();
    return rho_eff;
}

void viscosity_to_uniform_space(
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& b,
    const kynema_sgf::FieldRepo& repo,
    int lev)
{
    const auto& mesh_fac_xf =
        repo.get_mesh_mapping_field(kynema_sgf::FieldLoc::XFACE);
    const auto& mesh_fac_yf =
        repo.get_mesh_mapping_field(kynema_sgf::FieldLoc::YFACE);
    const auto& mesh_fac_zf =
        repo.get_mesh_mapping_field(kynema_sgf::FieldLoc::ZFACE);
    const auto& mesh_detJ_xf =
        repo.get_mesh_mapping_det_j(kynema_sgf::FieldLoc::XFACE);
    const auto& mesh_detJ_yf =
        repo.get_mesh_mapping_det_j(kynema_sgf::FieldLoc::YFACE);
    const auto& mesh_detJ_zf =
        repo.get_mesh_mapping_det_j(kynema_sgf::FieldLoc::ZFACE);

    // beta accounted for mesh mapping (x-face) = J/fac^2 * mu
    {
        const auto& mu_arrs = b[0].arrays();
        const auto& fac_arrs = mesh_fac_xf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_xf(lev).const_arrays();

        amrex::ParallelFor(
            b[0], [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                mu_arrs[nbx](i, j, k) =
                    mu_arrs[nbx](i, j, k) * detJ_arrs[nbx](i, j, k) /
                    kynema_sgf::utils::powi(fac_arrs[nbx](i, j, k, 0), 2);
            });
    }
    // beta accounted for mesh mapping (y-face) = J/fac^2 * mu
    {
        const auto& mu_arrs = b[1].arrays();
        const auto& fac_arrs = mesh_fac_yf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_yf(lev).const_arrays();

        amrex::ParallelFor(
            b[1], [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                mu_arrs[nbx](i, j, k) =
                    mu_arrs[nbx](i, j, k) * detJ_arrs[nbx](i, j, k) /
                    kynema_sgf::utils::powi(fac_arrs[nbx](i, j, k, 1), 2);
            });
    }
    // beta accounted for mesh mapping (z-face) = J/fac^2 * mu
    {
        const auto& mu_arrs = b[2].arrays();
        const auto& fac_arrs = mesh_fac_zf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_zf(lev).const_arrays();

        amrex::ParallelFor(
            b[2], [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                mu_arrs[nbx](i, j, k) =
                    mu_arrs[nbx](i, j, k) * detJ_arrs[nbx](i, j, k) /
                    kynema_sgf::utils::powi(fac_arrs[nbx](i, j, k, 2), 2);
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace diffusion
