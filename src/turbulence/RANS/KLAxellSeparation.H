#ifndef KLAXELLSEPARATION_H
#define KLAXELLSEPARATION_H

#include <memory>
#include <string>
#include "src/turbulence/RANS/KLAxell.H"

namespace kynema_sgf::turbulence {

namespace klaxell_separation {
//! Geometry component: height of the cell center above the surface
inline constexpr int geom_height = 0;
//! Geometry component: fluid weight of the cell (0 inside the terrain)
inline constexpr int geom_fluid_weight = 1;
//! Geometry component: Rt of the closure, after the neutral override
inline constexpr int geom_rt = 2;
//! Geometry component: curvature contribution to Rt (richardson model)
inline constexpr int geom_rt_curvature = 3;
//! Number of geometry components
inline constexpr int geom_ncomp = 4;
//! Dimensionless size below which the curvature correction does nothing,
//! so that round-off curvature in a straight flow leaves the model unchanged
inline constexpr amrex::Real curvature_tolerance = 1.0e-8_rt;
//! Deviation of the rotation function from 1 below which it does nothing; it
//! ramps in to the full deviation at twice this value. The strain-derivative
//! term is divided by D^4 and amplifies small departures from a horizontally
//! uniform flow where the shear is weak.
inline constexpr amrex::Real rotation_function_tolerance = 1.0e-3_rt;
} // namespace klaxell_separation

/** KLAxell one-equation RANS model with treatments for separating flow
 *  \ingroup turb_model
 *
 *  The closure is evaluated in a single kernel. The height above the
 *  surface and the fluid weight of each cell are filled beforehand by one
 *  geometry kernel per terrain treatment, so a new terrain treatment adds a
 *  geometry kernel instead of another copy of the closure. With every
 *  optional treatment disabled the model reproduces KLAxell.
 *
 *  Each treatment is switched on from the KLAxellSeparation input section.
 *  The pressure-gradient sensor is a diagnostic stored in the field
 *  pressure_gradient_sensor; the pressure-gradient treatments act only where
 *  the sensor exceeds its threshold. The curvature correction is not gated
 *  by the sensor and does nothing where the streamlines are straight.
 */
template <typename Transport>
class KLAxellSeparation : public KLAxell<Transport>
{
public:
    static std::string identifier()
    {
        return "KLAxellSeparation-" + Transport::identifier();
    }

    explicit KLAxellSeparation(CFDSim& sim);

    ~KLAxellSeparation() override;

    [[nodiscard]] std::string model_name() const override
    {
        return "KLAxellSeparation";
    }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! Actions to perform post initialization
    void post_init_actions() override;

    //! Actions to perform post regrid
    void post_regrid_actions() override;

    //! Relax the stored gate toward the latest sensor once per time step
    void post_advance_work() override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

    /** Fill the geometry fields for flat terrain
     *
     *  \param lev AMR level
     */
    void flat_geometry(int lev);

    /** Fill the geometry fields from the TerrainDrag fields
     *
     *  \param lev AMR level
     */
    void terrain_drag_geometry(int lev);

    /** Evaluate the closure on one level from the geometry fields
     *
     *  Also stores the Rt of the closure, after the neutral override, in the
     *  geometry field for the curvature correction.
     *
     *  \param lev AMR level
     *  \param fstate Field state
     *  \param beta Thermal expansion coefficient
     */
    void closure(int lev, FieldState fstate, const ScratchField& beta);

    /** Compute the pressure-gradient sensor on one level
     *
     *  The sensor is the pressure gradient along the local flow direction,
     *  (u / |u|) . grad(p), divided by rho (k + c_u |u|^2) / L. Inside the
     *  boundary layer the TKE sets the scale; above it, where the TKE
     *  vanishes, the velocity term keeps the sensor bounded.
     *
     *  \param lev AMR level
     *  \param fstate Field state
     */
    void pressure_gradient_sensor(int lev, FieldState fstate);

    /** Compute the pressure-gradient sensor on one level from the velocity
     *
     *  Along a streamline of steady inviscid flow the pressure gradient is
     *  -rho (u . grad)(|u|^2 / 2) / |u| (Bernoulli), so the sensor becomes
     *  -(u . grad)(|u|^2 / 2) / |u| divided by (k + c_u |u|^2) / L. The
     *  treatments leave the velocity smooth where they change the pressure
     *  field, so this estimate does not feed their changes back into the
     *  gate. The sensor is multiplied by the smallest fluid weight of the six
     *  face neighbors, which makes it zero in cells next to the terrain,
     *  where the velocity gradient reaches into the terrain.
     *
     *  \param lev AMR level
     *  \param fstate Field state
     */
    void velocity_sensor(int lev, FieldState fstate);

    /** Curvature correction through Rt on one level (richardson model)
     *
     *  The streamline curvature gives the Rayleigh / Bradshaw frequency
     *  N2_c = (2 / |u|^2) (|a_n|^2 - |u| a_n . grad|u|), with a = (u . grad) u
     *  and a_n its part normal to u; N2_c = 4 Omega^2 for solid rotation and
     *  0 for straight streamlines. Its contribution
     *  Rt_c = c_c L^2 N2_c / (Cmu^6 k) is added to the Rt of the closure, after
     *  its neutral override, and the eddy viscosity and the shear production
     *  are scaled by Cmu(Rt + Rt_c) / Cmu(Rt), the buoyancy production by
     *  Cmu'(Rt + Rt_c) / Cmu'(Rt). The length scale is unchanged. Rt_c is
     *  stored for the thermal diffusivity.
     *
     *  \param lev AMR level
     *  \param fstate Field state
     */
    void curvature_richardson(int lev, FieldState fstate);

    /** Curvature correction by the rotation function on one level
     *
     *  Spalart-Shur rotation function
     *  f = (1 + c_r1) (2 r* / (1 + r*)) (1 - c_r3 atan(c_r2 r~)) - c_r1,
     *  limited to [0, 1.25], with r* = S / Omega and
     *  r~ = 2 Omega_ik S_jk (DS_ij/Dt) / D^4, D^2 = (S^2 + Omega^2) / 2. The
     *  material derivative of the strain rate is taken as u . grad(S)
     *  (steady flow) and the frame rotation is neglected. The eddy viscosity
     *  and the shear and buoyancy production are multiplied by f; f = 1 in a
     *  simple shear and 0 in solid rotation.
     *
     *  \param lev AMR level
     *  \param fstate Field state
     */
    void curvature_rotation_function(int lev, FieldState fstate);

    /** Effective thermal diffusivity with the curvature contribution to Rt
     *
     *  Same as KLAxell::update_alphaeff, with the stored Rt_c of the
     *  richardson curvature model added to Rt for the turbulent Prandtl
     *  number.
     *
     *  \param alphaeff Effective thermal diffusivity
     */
    void update_alphaeff_curvature(Field& alphaeff);

    /** Apply the realizable Cmu limiter on one level
     *
     *  The eddy viscosity and the shear and buoyancy production are divided
     *  by 1 + c_s g max(0, Sigma / Cmu - 1), with Sigma = L S / sqrt(k) and
     *  the gate g = clamp((s - s_T) / s_T, 0, 1) of the pressure-gradient
     *  sensor s. The gate is exactly 0 below s_T, so those cells keep the
     *  KLAxell values. Sigma equals Cmu in an equilibrium log layer, and for
     *  large Sigma the eddy viscosity tends to rho Cmu(Rt) Cmu k / (c_s g S).
     *
     *  \param lev AMR level
     */
    void realizable_cmu(int lev);

    /** Apply the realizable Cmu limiter on one level with the stored gate
     *
     *  Same limiter as realizable_cmu, with the gate taken from the field
     *  separation_gate, which relax_gate updates once per time step.
     *
     *  \param lev AMR level
     */
    void realizable_cmu_relaxed(int lev);

    /** Relax the stored gate on one level over the last time step
     *
     *  The gate relaxes toward the ramp value g* = clamp((s - s_T) / s_T, 0, 1)
     *  of the latest sensor as g* + (g - g*) exp(-dt / tau), which is exact
     *  for a constant g* and stable for any time step. A gate at 0 stays
     *  exactly 0 while the sensor stays below s_T.
     *
     *  \param lev AMR level
     */
    void relax_gate(int lev);

protected:
    //! Height above the surface, fluid weight, Rt of the closure and
    //! curvature Rt of each cell, with one ghost cell for the neighbor weights
    std::unique_ptr<ScratchField> m_geometry;

    //! Compute the pressure-gradient sensor
    bool m_use_pressure_gradient_sensor{false};

    //! Weight c_u of the velocity term in the sensor scale
    amrex::Real m_sensor_velocity_weight{0.05_rt};

    //! Sensor value s_T where the treatments start; their weight reaches 1 at
    //! twice this value
    amrex::Real m_sensor_threshold{0.1_rt};

    //! Pressure-gradient sensor, declared only when the sensor is enabled
    Field* m_pressure_gradient_sensor{nullptr};

    //! Source of the sensor: "pressure" (the pressure gradient) or "velocity"
    //! (the Bernoulli estimate from the velocity field)
    std::string m_sensor_source{"pressure"};

    //! Velocity gradient, allocated only for the velocity sensor
    std::unique_ptr<ScratchField> m_grad_vel;

    //! Apply the realizable Cmu limiter
    bool m_use_realizable_cmu{false};

    //! Strength c_s of the realizable Cmu limiter
    amrex::Real m_realizable_cmu_strength{1.0_rt};

    //! Ratio C_P of the production cap applied by the TKE source
    amrex::Real m_production_cap_ratio{10.0_rt};

    //! Factor c_d on the dissipation of the destruction boost applied by the
    //! TKE source
    amrex::Real m_destruction_boost_factor{1.2_rt};

    //! Apply the curvature correction
    bool m_use_curvature_correction{false};

    //! Curvature model: "richardson" or "rotation_function"
    std::string m_curvature_model{"richardson"};

    //! Coefficient c_c on the curvature Rt of the richardson model
    amrex::Real m_curvature_coefficient{1.0_rt};

    //! Constants c_r1, c_r2 and c_r3 of the rotation function
    amrex::Real m_curvature_cr1{1.0_rt};
    amrex::Real m_curvature_cr2{12.0_rt};
    amrex::Real m_curvature_cr3{1.0_rt};

    //! Strain rate magnitude, kept for the realizable Cmu limiter
    std::unique_ptr<ScratchField> m_strain;

    //! Relaxation time tau of the gate (default 10 s); 0 uses the
    //! instantaneous gate
    amrex::Real m_gate_relaxation_time{10.0_rt};

    //! Gate relaxed in time, declared only when the sensor is on and the
    //! relaxation time is positive; kept in checkpoint files and
    //! interpolated on regrid
    Field* m_separation_gate{nullptr};
};

} // namespace kynema_sgf::turbulence

#endif /* KLAXELLSEPARATION_H */
