#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

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
    // amrex::Gpu::DeviceVector<amrex::Real> m_device_vel_ht;
    //   amrex::Gpu::DeviceVector<amrex::Real> m_device_vel_vals;
    std::string m_1d_rans;
    amrex::Vector<amrex::Real> m_wind_heights;
    amrex::Vector<amrex::Real> m_u_values;
    amrex::Vector<amrex::Real> m_v_values;
    amrex::Vector<amrex::Real> m_w_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_windht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_u_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_v_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_w_d;

    amrex::Real m_drag_coefficient{10.0};
    amrex::Real m_sponge_strength{1.0};
    amrex::Real m_sponge_density{1.0};
    amrex::Real m_sponge_distance_west{0};
    amrex::Real m_sponge_distance_east{0};
    amrex::Real m_sponge_distance_south{0};
    amrex::Real m_sponge_distance_north{0};
    int m_sponge_west{0};
    int m_sponge_east{0};
    int m_sponge_south{0};
    int m_sponge_north{0};
    bool m_is_laminar{false};
    bool m_terrain_is_waves{false};
    bool m_apply_MOSD{false};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
};

} // namespace amr_wind::pde::icns

#endif
