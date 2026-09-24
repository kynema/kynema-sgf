#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"
#include "src/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

#include <fstream>

using namespace amrex::literals;

namespace kynema_sgf::immersedterrain {

namespace {
//! Fractions within this tolerance of 0 or 1 are snapped, so that round-off in
//! the interpolated height still yields exactly fluid / solid cells and the
//! mask search has a clean threshold to work with.
constexpr amrex::Real fraction_tol = 1.0e-3_rt;
} // namespace

ImmersedTerrain::ImmersedTerrain(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_terrain_fraction(sim.repo().declare_field("terrain_fraction", 1, 1, 1))
    , m_terrain_mask(sim.repo().declare_int_field("terrain_mask", 1, 1, 1))
    , m_terrain_surface(sim.repo().declare_field("terrain_surface", 1, 1, 1))
{
    amrex::ParmParse pp(identifier());
    pp.query("terrain_file", m_terrain_file);
    pp.query("solid_threshold", m_solid_threshold);
    pp.query("drag_weight", m_drag_weight);
    if (m_drag_weight != "fraction" && m_drag_weight != "center") {
        amrex::Abort(
            identifier() + ".drag_weight must be fraction or center, got " +
            m_drag_weight);
    }

    {
        std::string turbulence_model{"Laminar"};
        amrex::ParmParse pp_turb("turbulence");
        pp_turb.query("model", turbulence_model);
        m_laminar = (turbulence_model == "Laminar");
    }

    pp.query("implicit_projection", m_implicit_projection);
    if (m_implicit_projection) {
        // Same coefficient the explicit source term would use
        amrex::ParmParse pp_drag("ImmersedDragForcing");
        pp_drag.query("drag_coefficient", m_drag_coefficient);
        m_terrain_drag_rate =
            &sim.repo().declare_field("terrain_drag_rate", 1, 1, 1);
        m_terrain_drag_rate->setVal(0.0_rt);
        m_terrain_drag_rate->set_default_fillpatch_bc(m_sim.time());
        m_sim.io_manager().register_io_var("terrain_drag_rate");
        amrex::Print() << identifier()
                       << ": immersed drag applied implicitly in the "
                          "projections with C_d = "
                       << m_drag_coefficient << "\n";
    }

    // No-slip wall flux at the terrain interface: face factors applied to
    // the diffusion coefficients of every equation (see pass 3)
    {
        const amrex::Array<FieldLoc, AMREX_SPACEDIM> locs{
            {FieldLoc::XFACE, FieldLoc::YFACE, FieldLoc::ZFACE}};
        const amrex::Array<std::string, AMREX_SPACEDIM> names{
            {"terrain_diffusion_xf", "terrain_diffusion_yf",
             "terrain_diffusion_zf"}};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            m_diffusion_factor[dir] =
                &sim.repo().declare_field(names[dir], 1, 0, 1, locs[dir]);
            m_diffusion_factor[dir]->setVal(1.0_rt);
        }
    }

    m_sim.io_manager().register_output_int_var("terrain_mask");
    m_sim.io_manager().register_io_var("terrain_fraction");
    m_sim.io_manager().register_io_var("terrain_surface");

    m_terrain_fraction.setVal(0.0_rt);
    m_terrain_mask.setVal(mask_fluid);
    m_terrain_surface.setVal(0.0_rt);
    m_terrain_fraction.set_default_fillpatch_bc(m_sim.time());
    m_terrain_surface.set_default_fillpatch_bc(m_sim.time());
}

void ImmersedTerrain::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::initialize_fields");

    //! Terrain coordinates from file
    amrex::Vector<amrex::Real> xterrain;
    amrex::Vector<amrex::Real> yterrain;
    amrex::Vector<amrex::Real> zterrain;
    ioutils::read_flat_grid_file(m_terrain_file, xterrain, yterrain, zterrain);

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& fraction = m_terrain_fraction(level);
    auto& mask = m_terrain_mask(level);
    auto& surface = m_terrain_surface(level);

    // Copy terrain to device
    const auto xterrain_size = xterrain.size();
    const auto yterrain_size = yterrain.size();
    const auto zterrain_size = zterrain.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(xterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_yterrain(yterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_zterrain(zterrain_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xterrain.begin(), xterrain.end(),
        d_xterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, yterrain.begin(), yterrain.end(),
        d_yterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, zterrain.begin(), zterrain.end(),
        d_zterrain.begin());
    const auto* xterrain_ptr = d_xterrain.data();
    const auto* yterrain_ptr = d_yterrain.data();
    const auto* zterrain_ptr = d_zterrain.data();

    auto frac_arrs = fraction.arrays();
    auto mask_arrs = mask.arrays();
    const bool has_rate = (m_terrain_drag_rate != nullptr);
    auto rate_arrs = has_rate ? (*m_terrain_drag_rate)(level).arrays()
                              : amrex::MultiArray4<amrex::Real>();
    const amrex::Real drag_rate_solid = m_drag_coefficient / dx[2];
    // Laminar flow with fraction weighting: the partial cell is a fluid cell
    // whose wall is carried by the no-slip wall flux, so only the cells
    // entirely inside the terrain are held by the drag
    const bool solid_only = m_laminar && (m_drag_weight == "fraction");
    const bool center_weight = (m_drag_weight == "center") || solid_only;
    const amrex::Real weight_threshold =
        solid_only ? 1.0_rt : m_solid_threshold;
    auto surf_arrs = surface.arrays();

    // Pass 1: terrain height and volume fraction, including
    // ghost cells so that the neighbor search in pass 2 has valid data.
    amrex::ParallelFor(
        fraction, m_terrain_fraction.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);

            const auto height = [=](const amrex::Real xq,
                                    const amrex::Real yq) {
                return interp::bilinear(
                    xterrain_ptr, xterrain_ptr + xterrain_size, yterrain_ptr,
                    yterrain_ptr + yterrain_size, zterrain_ptr, xq, yq);
            };

            const amrex::Real terrain_ht = height(x, y);
            surf_arrs[nbx](i, j, k, surf_height) = terrain_ht;

            // Fraction of the cell column below the terrain height
            const amrex::Real z_bottom = prob_lo[2] + (k * dx[2]);
            amrex::Real vol_frac = (terrain_ht - z_bottom) / dx[2];
            vol_frac = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(vol_frac, 0.0_rt), 1.0_rt);
            if (vol_frac < fraction_tol) {
                vol_frac = 0.0_rt;
            } else if (vol_frac > 1.0_rt - fraction_tol) {
                vol_frac = 1.0_rt;
            }
            // Ghost cells below the domain floor stay fluid so the bottom row
            // is not flagged as a terrain surface by the neighbor search
            frac_arrs[nbx](i, j, k, 0) = (z > prob_lo[2]) ? vol_frac : 0.0_rt;
            if (has_rate) {
                rate_arrs[nbx](i, j, k, 0) =
                    kynema_sgf::immersed_wall::solid_weight(
                        frac_arrs[nbx](i, j, k, 0), center_weight,
                        weight_threshold) *
                    drag_rate_solid;
            }
        });
    amrex::Gpu::streamSynchronize();
    // Ghost cells across periodic boundaries take the wrapped values instead
    // of the clamped interpolation of the terrain file; physical-boundary
    // ghosts keep the analytic values from pass 1
    fraction.FillBoundary(geom.periodicity());
    surface.FillBoundary(geom.periodicity());
    if (has_rate) {
        (*m_terrain_drag_rate)(level).FillBoundary(geom.periodicity());
    }

    // Pass 2: cell classification. A fluid cell becomes a surface cell if any
    // of its six face neighbors is a wall (see immersed_wall::wall_threshold);
    // this catches the side walls of steep terrain and buildings that a
    // vertical-only search misses.
    const amrex::Real wall_fraction = immersed_wall::wall_threshold(
        m_drag_weight == "center", m_solid_threshold);
    amrex::ParallelFor(
        fraction, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& frac = frac_arrs[nbx];
            const amrex::Real f = frac(i, j, k, 0);
            int cell_mask = mask_fluid;
            if (f >= 1.0_rt) {
                cell_mask = mask_solid;
            } else if (f > 0.0_rt) {
                cell_mask = mask_surface;
            } else {
                const amrex::Real max_nb = amrex::max<amrex::Real>(
                    amrex::max<amrex::Real>(
                        frac(i - 1, j, k, 0), frac(i + 1, j, k, 0)),
                    amrex::max<amrex::Real>(
                        frac(i, j - 1, k, 0), frac(i, j + 1, k, 0)),
                    amrex::max<amrex::Real>(
                        frac(i, j, k - 1, 0), frac(i, j, k + 1, 0)));
                if (max_nb >= wall_fraction) {
                    cell_mask = mask_surface;
                }
            }
            mask_arrs[nbx](i, j, k, 0) = cell_mask;
        });
    amrex::Gpu::streamSynchronize();

    // Pass 3: face factors for the diffusion coefficients. On a wall face the
    // coefficient is scaled by dx_f / d1, so that the discrete flux
    // mu (u - 0) / dx_f equals the no-slip flux mu u / d1 at the true wall
    // position. With fraction weighting the partial cell stands for the
    // centroid of its fluid part: its bottom face uses d1 = (1 - beta) dz / 2
    // and the face above it spans (2 - beta) dz / 2. The fraction terms vanish
    // with center weighting, where d1 = z - h.
    const amrex::Real centroid =
        (m_drag_weight == "fraction") ? 1.0_rt : 0.0_rt;
    // The terrain height is resolved to fraction_tol * dz (pass 1), so d1 is
    // the true wall distance down to that resolution; the floor only keeps
    // the factor finite when the wall falls on a cell center
    const amrex::Real d1_min = 0.5_rt * fraction_tol;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        auto& ff = (*m_diffusion_factor[dir])(level);
        auto fac_arrs = ff.arrays();
        amrex::ParallelFor(
            ff, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // Face (i,j,k) in direction dir separates the low cell
                // (i,j,k) - e_dir from the high cell (i,j,k)
                const int il = i - ((dir == 0) ? 1 : 0);
                const int jl = j - ((dir == 1) ? 1 : 0);
                const int kl = k - ((dir == 2) ? 1 : 0);
                const auto& frac = frac_arrs[nbx];
                const amrex::Real beta_lo = frac(il, jl, kl, 0);
                const amrex::Real beta_hi = frac(i, j, k, 0);
                const bool solid_lo = beta_lo >= wall_fraction;
                const bool solid_hi = beta_hi >= wall_fraction;
                amrex::Real factor = 1.0_rt;
                if (solid_lo != solid_hi) {
                    // Wall at the face unless the terrain height gives the
                    // true distance below a fluid cell
                    amrex::Real d1 = 0.5_rt * dx[dir];
                    if (dir == 2 && solid_lo) {
                        const amrex::Real z_c =
                            prob_lo[2] + ((k + 0.5_rt) * dx[2]);
                        const amrex::Real h = surf_arrs[nbx](
                            i, j, k, ImmersedTerrain::surf_height);
                        const amrex::Real shift =
                            centroid * (beta_hi > 0.0_rt ? 1.0_rt : 0.0_rt);
                        d1 = (z_c - h) +
                             (shift * 0.5_rt * ((0.5_rt * dx[2]) - (z_c - h)));
                        d1 = amrex::min<amrex::Real>(
                            amrex::max<amrex::Real>(d1, d1_min * dx[2]), dx[2]);
                    }
                    factor = dx[dir] / d1;
                } else if (dir == 2 && !solid_lo) {
                    // Fluid face above a partial cell
                    factor = 2.0_rt / (2.0_rt - (centroid * beta_lo));
                }
                fac_arrs[nbx](i, j, k, 0) = factor;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void ImmersedTerrain::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

} // namespace kynema_sgf::immersedterrain
