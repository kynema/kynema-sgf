#ifndef VORTEXRING_H
#define VORTEXRING_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/** Vortex ring physics
 *  \ingroup physics
 */
class VortexRing : public Physics::Register<VortexRing>
{
public:
    static std::string identifier() { return "VortexRing"; }

    explicit VortexRing(const CFDSim& sim);

    ~VortexRing() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

    template <typename VortexRingType>
    void initialize_velocity(const VortexRingType& vorticity_theta);

private:
    const CFDSim& m_sim;
    const FieldRepo& m_repo;
    Field& m_velocity;
    Field& m_density;

    //! initial density value
    amrex::Real m_rho{1.0_rt};

    amrex::Real m_R{1.0_rt};
    amrex::Real m_Gamma{1.0_rt};
    amrex::Real m_delta{0.2_rt};
    amrex::Real m_dz{2.0_rt};
    amrex::Vector<int> m_perturbation_modes;
    amrex::Vector<amrex::Real> m_perturbation_phases_1;
    amrex::Vector<amrex::Real> m_perturbation_phases_2;
    // amplitude relative to the unperturbed vortex ring diameter
    amrex::Real m_perturbation_amplitude{1.0e-4_rt};

    std::string m_vortexringtype{"fatcore"};
};

struct FatCore
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*r*/,
        amrex::Real /*unused*/,
        amrex::Real /*z*/,
        amrex::Real /*R*/,
        amrex::Real /*Gamma*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        amrex::Real /*unused*/,
        int /*unused*/,
        const int* /*unused*/,
        const amrex::Real* /*unused*/,
        const amrex::Real* /*unused*/) const;
};

struct CollidingRings
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real operator()(
        amrex::Real /*r*/,
        amrex::Real /*theta*/,
        amrex::Real /*z*/,
        amrex::Real /*R*/,
        amrex::Real /*Gamma*/,
        amrex::Real /*delta*/,
        amrex::Real /*dz*/,
        amrex::Real /*perturbation_amplitude*/,
        int /*num_modes*/,
        const int* /*perturbation_modes*/,
        const amrex::Real* /*perturbation_phases_1*/,
        const amrex::Real* /*perturbation_phases_2*/) const;
};

} // namespace amr_wind

#endif /* VortexRing_H */
