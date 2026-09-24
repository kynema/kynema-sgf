#include "src/equation_systems/icns/source_terms/ImmersedDragForcing.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
using kynema_sgf::immersedterrain::ImmersedTerrain;

namespace kynema_sgf::pde::icns {

ImmersedDragForcing::ImmersedDragForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp(identifier());
    pp.query("drag_coefficient", m_drag_coefficient);

    amrex::ParmParse pp_terrain(ImmersedTerrain::identifier());
    pp_terrain.query("solid_threshold", m_solid_threshold);
    pp_terrain.query("drag_weight", m_drag_weight);

    std::string turbulence_model = "Laminar";
    amrex::ParmParse pp_turb("turbulence");
    pp_turb.query("model", turbulence_model);
    m_is_laminar = (turbulence_model == "Laminar");

    if (!sim.repo().field_exists("terrain_fraction")) {
        amrex::Abort(
            identifier() +
            " requires the ImmersedTerrain physics "
            "(terrain_fraction field not found)");
    }
}

ImmersedDragForcing::~ImmersedDragForcing() = default;

void ImmersedDragForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    const auto& repo = m_sim.repo();
    // With ImmersedTerrain.implicit_projection the drag is applied inside the
    // projections through terrain_drag_rate
    if (repo.field_exists("terrain_drag_rate")) {
        return;
    }

    auto const& src_arrs = src_term.arrays();
    // C_eff integrates the drag exactly over the whole step from u^n, so the
    // source always acts on the old velocity: evaluated on the predictor
    // velocity in a corrector it would give 1 - (1 - e) e instead of
    // e = exp(-C dt)
    auto const& vel_arrs =
        m_velocity.state(FieldState::Old)(lev).const_arrays();
    auto const& frac_arrs =
        repo.get_field("terrain_fraction")(lev).const_arrays();

    const auto& geom = m_mesh.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const amrex::Real dt = m_time.delta_t();
    const amrex::Real drag_rate = m_drag_coefficient / dx[2];

    // Laminar with fraction weighting: the partial cell is a fluid cell and
    // its wall comes from the no-slip wall flux, so only cells entirely
    // inside the terrain take the drag
    const bool center_weight = (m_drag_weight == "center");
    const bool solid_only = m_is_laminar && !center_weight;
    const bool drag_center = center_weight || solid_only;
    const amrex::Real drag_threshold = solid_only ? 1.0_rt : m_solid_threshold;

    amrex::ParallelFor(
        src_term, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real w_solid = kynema_sgf::immersed_wall::solid_weight(
                frac_arrs[nbx](i, j, k, 0), drag_center, drag_threshold);
            if (w_solid <= 0.0_rt) {
                return;
            }
            const amrex::Real C_eff =
                kynema_sgf::immersed_wall::exact_relaxation_rate(
                    w_solid * drag_rate, dt);
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                src_arrs[nbx](i, j, k, n) -= C_eff * vel_arrs[nbx](i, j, k, n);
            }
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace kynema_sgf::pde::icns
