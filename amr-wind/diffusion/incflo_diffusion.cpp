#include "amr-wind/incflo.H"
#include "amr-wind/diffusion/diffusion.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace diffusion {

amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
get_diffuse_tensor_bc(amr_wind::Field& velocity, amrex::Orientation::Side side)
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
get_diffuse_scalar_bc(amr_wind::Field& scalar, amrex::Orientation::Side side)
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
    BL_PROFILE("amr-wind::diffusion::average_velocity_eta_to_faces");
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
    BL_PROFILE("amr-wind::diffusion::fixup_eta_on_domain_faces");

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

void viscosity_to_uniform_space(
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& b,
    const amr_wind::FieldRepo& repo,
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

    // beta accounted for mesh mapping (x-face) = J/fac^2 * mu
    {
        const auto& mu_arrs = b[0].arrays();
        const auto& fac_arrs = mesh_fac_xf(lev).arrays();
        const auto& detJ_arrs = mesh_detJ_xf(lev).const_arrays();

        amrex::ParallelFor(
            b[0], [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                mu_arrs[nbx](i, j, k) =
                    mu_arrs[nbx](i, j, k) * detJ_arrs[nbx](i, j, k) /
                    std::pow(fac_arrs[nbx](i, j, k, 0), 2.0_rt);
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
                    std::pow(fac_arrs[nbx](i, j, k, 1), 2.0_rt);
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
                    std::pow(fac_arrs[nbx](i, j, k, 2), 2.0_rt);
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace diffusion
