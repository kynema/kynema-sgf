#ifndef ABLMEANBOUSSINESQ_H
#define ABLMEANBOUSSINESQ_H

#include "src/core/FieldRepo.H"
#include "src/transport_models/TransportModel.H"
#include "src/equation_systems/icns/MomentumSource.H"
#include "src/utilities/FieldPlaneAveraging.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

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
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

    void mean_temperature_init(const FieldPlaneAveraging& /*tavg*/);

    void mean_temperature_update(const FieldPlaneAveraging& /*tavg*/);

private:
    const amrex::AmrCore& m_mesh;
    FieldRepo& m_repo;

    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    amrex::Gpu::DeviceVector<amrex::Real> m_theta_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_vals;

    //! Transport model
    const transport::TransportModel& m_transport;

    int m_axis{2};

    bool m_const_profile{false};

    mutable std::unique_ptr<ScratchField> m_beta_scratch;
    mutable std::unique_ptr<ScratchField> m_ref_theta_scratch;

    //! Read a temperature profile from a text file
    void read_temperature_profile(const std::string& profile_file_name);
};
} // namespace kynema_sgf::pde::icns

#endif /* ABLMEANBOUSSINESQ_H */
