#ifndef KLAXELL_H
#define KLAXELL_H

#include <string>
#include "src/core/ScratchField.H"
#include "src/turbulence/TurbModelBase.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::turbulence {

/** Base class for 1-Equation RANS TKE turbulence model
 *  \ingroup turb_model
 */
template <typename Transport>
class KLAxell : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "KLAxell-" + Transport::identifier();
    }

    explicit KLAxell(CFDSim& sim);

    [[nodiscard]] std::string model_name() const override { return "KLAxell"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! Do any post advance work
    void post_advance_work() override;

    //! Actions to perform post initialization
    void post_init_actions() override;

    //! Actions to perform post regrid
    void post_regrid_actions() override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Update the effective scalar diffusivity field
    void update_scalar_diff(Field& deff, const std::string& name) override;

    /** With TerrainDrag and ``KLAxell.terrain_wall_stencil``: strain rate of
     *  the fluid cells next to the blanked cells from the stencil that flat
     *  ground uses at its bottom wall, the wall being the face of the blanked
     *  cell. Public because it launches a device lambda.
     *
     *  @param velocity Velocity field of the state being updated
     */
    void terrain_wall_strain_rate(const Field& velocity);

    /** With TerrainDrag and ``KLAxell.terrain_blanked_face_length``: height
     *  of the top face of the blanked column, from which the mixing length is
     *  measured in place of the terrain height. Public because it launches a
     *  device lambda.
     *
     *  @param face_height Scratch field (one component) set to the face height
     */
    void terrain_blanked_face_height(ScratchField& face_height) const;

    /** With TerrainDrag and ``KLAxell.terrain_face_stress``: viscosity of the
     *  drag cells sized so that the face above each one carries the wall
     *  stress u*^2 of the DragForcing wall law. Public because it launches a
     *  device lambda.
     *
     *  @param lev Level to update
     *  @param velocity Velocity field of the state being updated
     *  @param density Density field of the state being updated
     */
    void
    terrain_face_stress(int lev, const Field& velocity, const Field& density);

    /** With TerrainDrag, ``KLAxell.terrain_face_heat_flux`` and
     *  ``ABL.wall_het_model = mol``: heat diffusivity of the drag cells sized
     *  so that the face above each one carries the surface heat flux
     *  -u* theta* given by the Obukhov length. Public because it launches a
     *  device lambda.
     *
     *  @param lev Level to update
     *  @param alphaeff Effective thermal diffusivity being updated
     *  @param lam_alpha Laminar thermal diffusivity
     */
    void terrain_face_heat_flux(
        int lev, Field& alphaeff, const ScratchField& lam_alpha);

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    Field& m_vel;
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_rho;
    Field* m_tke{nullptr};
    //! Turbulence constant
    amrex::Real m_Cmu{0.556_rt};
    amrex::Real m_Cmu_prime{0.556_rt};
    amrex::Real m_Cb_stable{0.25_rt};
    amrex::Real m_Cb_unstable{0.35_rt};
    amrex::Real m_prandtl{1.0_rt};
    Field& m_temperature;
    //! Gravity vector (m/s^2)
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
    amrex::Real m_surf_flux{0};
    amrex::Real m_meso_sponge_start{2000};

    std::unique_ptr<ScratchField> m_gradT;

    //! TerrainDrag options (``KLAxell.<name>``), all off by default and all
    //! on with ``TerrainDrag.wall_treatment = improved``
    //! Flat-ground wall stencil for the strain rate of the wall cells
    bool m_terrain_wall_stencil{false};
    //! Mixing length measured from the top face of the blanked column
    bool m_terrain_blanked_face_length{false};
    //! Drag-cell viscosity from the wall stress on the face above it
    bool m_terrain_face_stress{false};
    //! Drag-cell heat diffusivity from the surface heat flux on the face
    //! above it (mol only)
    bool m_terrain_face_heat_flux{false};

    //! Surface-layer inputs of the face stress and heat flux, read from the
    //! same inputs as DragForcing (``ABL.wall_het_model``,
    //! ``ABL.monin_obukhov_length``, ``ABL.kappa``, ``ABL.mo_beta_m``,
    //! ``ABL.mo_gamma_m``, ``DragForcing.minimum_z0``)
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_beta_m{16.0_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_min_z0{1.0e-4_rt};
};

} // namespace kynema_sgf::turbulence

#endif /* KLAXELL_H */
