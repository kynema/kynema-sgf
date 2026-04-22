#include "amr-wind/overset/overset_ops_routines.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::overset_ops {

/** Convert iblanks to AMReX mask
 *
 *  \f{align}
 *  \mathrm{mask}_{i,j,k} = \begin{cases}
 *  1 & \mathrm{IBLANK}_{i, j, k} = 1 \\
 *  0 & \mathrm{IBLANK}_{i, j, k} \leq 0
 *  \end{cases}
 *  \f}
 */
void iblank_to_mask(const IntField& iblank, IntField& maskf)
{
    const auto& nlevels = iblank.repo().mesh().finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ibl = iblank(lev);
        auto& mask = maskf(lev);

        const auto& ibarrs = ibl.const_arrays();
        const auto& marrs = mask.arrays();
        amrex::ParallelFor(
            ibl, ibl.n_grow,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                marrs[nbx](i, j, k) = amrex::max(ibarrs[nbx](i, j, k), 0);
            });
    }
    amrex::Gpu::streamSynchronize();
}

/** VOF-sensitive conversion of iblank to mask at nodes
 *
 * Masks most of the domain, including a layer around solid bodies, but avoids
 * masking cells near the interface
 */
void iblank_node_to_mask_vof(
    const IntField& iblank, const Field& voff, IntField& maskf)
{
    const auto& nlevels = iblank.repo().mesh().finestLevel() + 1;
    constexpr amrex::Real band_tol = 1.0e-4_rt;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ibl = iblank(lev);
        const auto& vof = voff(lev);
        auto& mask = maskf(lev);

        const auto& ibarrs = ibl.const_arrays();
        const auto& vofarrs = vof.const_arrays();
        const auto& marrs = mask.arrays();
        amrex::ParallelFor(
            ibl, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // Default is masking all 0 and -1 iblanks
                marrs[nbx](i, j, k) = amrex::max(ibarrs[nbx](i, j, k), 0);
                // Check cells neighboring node for being near interface
                bool near_interface = false;
                for (int ii = i - 1; ii < i + 1; ii++) {
                    for (int jj = j - 1; jj < j + 1; jj++) {
                        for (int kk = k - 1; kk < k + 1; kk++) {
                            near_interface =
                                near_interface ||
                                amr_wind::multiphase::interface_band(
                                    ii, jj, kk, vofarrs[nbx], 1, band_tol);
                        }
                    }
                }
                // Check nodes neighboring node for being near solid or overset
                bool near_solid = false;
                bool near_overset_boundary = false;
                for (int ii = i - 1; ii < i + 2; ii++) {
                    for (int jj = j - 1; jj < j + 2; jj++) {
                        for (int kk = k - 1; kk < k + 2; kk++) {
                            near_solid = ibarrs[nbx](ii, jj, kk) == 0
                                             ? true
                                             : near_solid;
                            near_overset_boundary = ibarrs[nbx](ii, jj, kk) == 1
                                                        ? true
                                                        : near_overset_boundary;
                        }
                    }
                }
                // Do mask -1 cells near interface and overset boundary
                if (ibarrs[nbx](i, j, k) == -1 &&
                    (near_overset_boundary && near_interface && !near_solid)) {
                    marrs[nbx](i, j, k) = 1;
                }
            });
        mask.FillBoundary(maskf.repo().mesh().Geom(lev).periodicity());
    }
    amrex::Gpu::streamSynchronize();
}

/** VOF-sensitive conversion of iblank to mask at cells
 *
 * Same concept as iblank_node_to_mask_vof, but for MAC projection
 */
void prepare_mask_cell_for_mac(FieldRepo& repo)
{
    const bool vof_exists = repo.field_exists("vof");

    if (vof_exists) {
        IntField& mask_field = repo.get_int_field("mask_cell");
        const IntField& iblank = repo.get_int_field("iblank_cell");
        const Field& f_vof = repo.get_field("vof");
        const auto& nlevels = repo.mesh().finestLevel() + 1;
        constexpr amrex::Real band_tol = 1.0e-4_rt;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& ibl = iblank(lev);
            const auto& vof = f_vof(lev);
            auto& mask = mask_field(lev);

            const auto& ibarrs = ibl.const_arrays();
            const auto& vofarrs = vof.const_arrays();
            const auto& marrs = mask.arrays();
            amrex::ParallelFor(
                ibl, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    // Default is masking all 0 and -1 iblanks
                    marrs[nbx](i, j, k) = amrex::max(ibarrs[nbx](i, j, k), 0);
                    // Check cells neighboring node for being near interface
                    bool near_interface = amr_wind::multiphase::interface_band(
                        i, j, k, vofarrs[nbx], 1, band_tol);
                    // Check neighboring cells for being near solid or overset
                    bool near_solid = false;
                    bool near_overset_boundary = false;
                    for (int ii = i - 1; ii < i + 2; ii++) {
                        for (int jj = j - 1; jj < j + 2; jj++) {
                            for (int kk = k - 1; kk < k + 2; kk++) {
                                near_solid = ibarrs[nbx](ii, jj, kk) == 0
                                                 ? true
                                                 : near_solid;
                                near_overset_boundary =
                                    ibarrs[nbx](ii, jj, kk) == 1
                                        ? true
                                        : near_overset_boundary;
                            }
                        }
                    }
                    // Do mask -1 cells near interface
                    if (ibarrs[nbx](i, j, k) == -1 &&
                        (near_overset_boundary && near_interface &&
                         !near_solid)) {
                        marrs[nbx](i, j, k) = 1;
                    }
                });
            mask.FillBoundary(mask_field.repo().mesh().Geom(lev).periodicity());
        }
        amrex::Gpu::streamSynchronize();
    }
}

/** Convert iblank to mask using ordinary method
 *
 * Intended to return mask_cell back to original values following
 * prepare_mask_cell_for_mac because mask_cell is used in other parts of the
 * flow solver, not just the MAC projection
 */
void revert_mask_cell_after_mac(FieldRepo& repo)
{
    const bool vof_exists = repo.field_exists("vof");

    if (vof_exists) {
        IntField& maskf = repo.get_int_field("mask_cell");
        const IntField& iblank = repo.get_int_field("iblank_cell");

        iblank_to_mask(iblank, maskf);
    }
}

// Swap pressure gradient values in overset region
void replace_gradp(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_gp0,
    const amrex::iMultiFab& mf_iblank)
{
    const auto gp = mf_gp.arrays();
    const auto gp0 = mf_gp0.const_arrays();
    const auto iblank = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_gp, mf_gp.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            if (iblank[nbx](i, j, k) <= 0) {
                gp[nbx](i, j, k, 0) = gp0[nbx](i, j, k, 0);
                gp[nbx](i, j, k, 1) = gp0[nbx](i, j, k, 1);
                gp[nbx](i, j, k, 2) = gp0[nbx](i, j, k, 2);
            }
        });
}

// Apply pressure gradient to velocity field
void apply_pressure_gradient(
    amrex::MultiFab& mf_vel,
    const amrex::MultiFab& mf_density,
    const amrex::MultiFab& mf_gp,
    const amrex::Real scaling_factor)
{
    const auto& vel = mf_vel.arrays();
    const auto& rho = mf_density.const_arrays();
    const auto& gp = mf_gp.const_arrays();
    amrex::ParallelFor(
        mf_vel, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real soverrho = scaling_factor / rho[nbx](i, j, k);
            vel[nbx](i, j, k, 0) -= gp[nbx](i, j, k, 0) * soverrho;
            vel[nbx](i, j, k, 1) -= gp[nbx](i, j, k, 1) * soverrho;
            vel[nbx](i, j, k, 2) -= gp[nbx](i, j, k, 2) * soverrho;
        });
}

} // namespace amr_wind::overset_ops
