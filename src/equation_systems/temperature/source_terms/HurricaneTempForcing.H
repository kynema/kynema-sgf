#ifndef HURRICANE_TEMP_FORCING_H
#define HURRICANE_TEMP_FORCING_H

#include "src/core/FieldRepo.H"
#include "src/equation_systems/temperature/TemperatureSource.H"
#include "src/utilities/FieldPlaneAveraging.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

class HurricaneTempForcing
    : public TemperatureSource::Register<HurricaneTempForcing>
{

public:
    static std::string identifier() { return "HurricaneTempForcing"; }

    explicit HurricaneTempForcing(const CFDSim& /*sim*/);

    ~HurricaneTempForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& /*mfi*/,
        const amrex::Box& bx,
        FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void mean_velocity_init(const VelPlaneAveraging& /*vavg*/);

    void mean_velocity_update(const VelPlaneAveraging& /*vavg*/);

private:
    //! Mesh
    const amrex::AmrCore& m_mesh;

    amrex::Gpu::DeviceVector<amrex::Real> m_vel_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_vel_vals;

    //! Axis over which averages are computed
    int m_axis{2};

    //! Temperature radial decay
    amrex::Real m_dTdR{0.001_rt};

    //! Zero Temperature height
    amrex::Real m_dTzh{18000.};
};

} // namespace kynema_sgf::pde::temperature

#endif
