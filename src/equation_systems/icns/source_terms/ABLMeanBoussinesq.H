#ifndef ABLMEANBOUSSINESQ_H
#define ABLMEANBOUSSINESQ_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Boussinesq buoyancy source term
 *  \ingroup icns_src we_abl
 *
 *  \f[
 *    S = \beta g \left( T_\mathrm{mean} - T_\mathrm{ref} \right)
 *  \f]
 */
class ABLMeanBoussinesq : public MomentumSource::Register<ABLMeanBoussinesq>
{
public:
    static std::string identifier() { return "ABLMeanBoussinesq"; }

    explicit ABLMeanBoussinesq(const CFDSim& sim);

    ~ABLMeanBoussinesq() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void mean_temperature_init(const FieldPlaneAveraging& /*tavg*/);

    void mean_temperature_update(const FieldPlaneAveraging& /*tavg*/);

private:
    const amrex::AmrCore& m_mesh;

    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    amrex::Gpu::DeviceVector<amrex::Real> m_theta_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_vals;

    //! Transport model
    const transport::TransportModel& m_transport;

    int m_axis{2};

    bool m_const_profile{false};

    //! Read a temperature profile from a text file
    void read_temperature_profile(const std::string& profile_file_name);
};
} // namespace amr_wind::pde::icns

#endif /* ABLMEANBOUSSINESQ_H */
