#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/constants.H"

namespace amr_wind::pde::icns {

/** Adds the forcing term to include the presence of immersed boundary
 *
 *  \ingroup icns_src
 *
 *
 */
class DragForcing : public MomentumSource::Register<DragForcing>
{
public:
    static std::string identifier() { return "DragForcing"; }

    explicit DragForcing(const CFDSim& sim);

    ~DragForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field* m_target_vel{nullptr};
    const Field* m_target_levelset{nullptr};
    amrex::Real m_drag_coefficient{10.0};
    bool m_is_laminar{false};
    bool m_terrain_is_waves{false};
    bool m_apply_MOSD{false};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41};
    amrex::Real m_gamma_m{5.0};
    amrex::Real m_beta_m{16.0};
};

} // namespace amr_wind::pde::icns

#endif
