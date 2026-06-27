#ifndef PARTIALDRAGFORCING_H
#define PARTIALDRAGFORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Adds the forcing term to include the presence of partial terrain drag/wall
 *
 *  \ingroup icns_src
 *
 *
 */
class PartialDragForcing : public MomentumSource::Register<PartialDragForcing>
{
public:
    static std::string identifier() { return "PartialDragForcing"; }

    explicit PartialDragForcing(const CFDSim& sim);

    ~PartialDragForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    amrex::Real m_bc_forcing_time_factor{5.0_rt};
    amrex::Real m_drag_coefficient{10.0_rt};
    amrex::Real m_min_z{0.1_rt};
    int m_drag_forcing{0};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
};

} // namespace kynema_sgf::pde::icns

#endif
