#ifndef RAYLEIGH_DAMPING_H
#define RAYLEIGH_DAMPING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Apply Rayleigh damping as a source terms
 * \ingroup icns_src
 *
 * The implementation is based on a private communication I (G. Deskos)
 * had with Matt Churchfield
 */

class RayleighDamping : public MomentumSource::Register<RayleighDamping>
{
public:
    static std::string identifier() { return "RayleighDamping"; }

    explicit RayleighDamping(const CFDSim& /*sim*/);

    ~RayleighDamping() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const amrex::AmrCore& m_mesh;

    const Field& m_velocity;

    //! Reference velocity defined as an input
    amrex::Vector<amrex::Real> m_ref_vel{15.0_rt, 0.0_rt, 0.0_rt};

    //! Full length of the damping layer
    amrex::Real m_dRD{500.};

    //! Damping length
    amrex::Real m_dFull{300.};

    //! Time scale defined by the user
    amrex::Real m_tau{5.0_rt};

    //! Which coordinate directions are forced
    amrex::Vector<int> m_fcoord{1, 1, 1};
};

} // namespace amr_wind::pde::icns

#endif
