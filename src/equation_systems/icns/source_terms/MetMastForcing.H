#ifndef METMASTFORCING_H
#define METMASTFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

class MetMastForcing : public MomentumSource::Register<MetMastForcing>
{
public:
    static std::string identifier() { return "MetMastForcing"; }

    explicit MetMastForcing(const CFDSim& sim);

    ~MetMastForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    std::string m_1d_metmast;
    const Field& m_velocity;
    amrex::Vector<amrex::Real> m_metmast_x;
    amrex::Vector<amrex::Real> m_metmast_y;
    amrex::Vector<amrex::Real> m_metmast_z;
    amrex::Vector<amrex::Real> m_u_values;
    amrex::Vector<amrex::Real> m_v_values;
    amrex::Vector<amrex::Real> m_w_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_x_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_y_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_z_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_u_values_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_v_values_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_w_values_d;
    amrex::Real m_meso_timescale{30};
    amrex::Real m_long_radius{500};
    amrex::Real m_vertical_radius{25};
    amrex::Real m_damping_radius{1400};
    const CFDSim& m_sim;
};

} // namespace amr_wind::pde::icns
#endif
