#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

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
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field* m_target_vel{nullptr};
    const Field* m_target_levelset{nullptr};
    std::string m_1d_rans;
    amrex::Vector<amrex::Real> m_wind_heights;
    amrex::Vector<amrex::Real> m_u_values;
    amrex::Vector<amrex::Real> m_v_values;
    amrex::Vector<amrex::Real> m_w_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_windht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_u_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_v_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_w_d;

    amrex::Real m_drag_coefficient{10.0_rt};
    amrex::Real m_sponge_strength{1.0_rt};
    amrex::Real m_sponge_density{1.0_rt};
    amrex::Real m_forcing_time_factor{5.0_rt};
    amrex::Real m_sponge_distance_west{0.0_rt};
    amrex::Real m_sponge_distance_east{0.0_rt};
    amrex::Real m_sponge_distance_south{0.0_rt};
    amrex::Real m_sponge_distance_north{0.0_rt};
    amrex::Real m_min_z{0.1_rt};
    bool m_sponge_west{false};
    bool m_sponge_east{false};
    bool m_sponge_south{false};
    bool m_sponge_north{false};
    bool m_is_laminar{false};
    bool m_terrain_is_waves{false};
    bool m_apply_MOSD{false};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
};

} // namespace kynema_sgf::pde::icns

#endif
