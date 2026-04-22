#ifndef CORIOLISFORCING_H
#define CORIOLISFORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Account for Coriolis forcing due to earth's rotation
 *  \ingroup icns_src we_abl
 *
 */
class CoriolisForcing : public MomentumSource::Register<CoriolisForcing>
{
public:
    static std::string identifier() { return "CoriolisForcing"; }

    explicit CoriolisForcing(const CFDSim& /*sim*/);

    ~CoriolisForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_velocity;

    ///@{
    /** Orientation of the grid w.r.t. to planetary coordinate system
     */
    amrex::Vector<amrex::Real> m_east{1.0_rt, 0.0_rt, 0.0_rt};
    amrex::Vector<amrex::Real> m_north{0.0_rt, 1.0_rt, 0.0_rt};
    amrex::Vector<amrex::Real> m_up{0.0_rt, 0.0_rt, 1.0_rt};
    ///@}

    //! Latitude where the Coriolis forcing is computed (+ve northern
    //! hemisphere)
    amrex::Real m_latitude{0.0_rt};

    amrex::Real m_sinphi{0.0_rt};
    amrex::Real m_cosphi{0.0_rt};

    //! omega --> Earth’s rotation angular speed
    amrex::Real m_omega{7.3e-5_rt};

    //! 2.0_rt * omega
    amrex::Real m_coriolis_factor{0.0_rt};

    bool m_is_horizontal{false};
};

} // namespace amr_wind::pde::icns

#endif /* CORIOLISFORCING_H */
