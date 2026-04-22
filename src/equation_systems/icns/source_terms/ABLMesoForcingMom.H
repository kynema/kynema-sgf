#ifndef ABLMESOFORCINGMOM_H
#define ABLMESOFORCINGMOM_H

#include <memory>
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/wind_energy/ABLMesoscaleForcing.H"
#include "amr-wind/wind_energy/ABLMesoscaleInput.H"

namespace amr_wind::pde::icns {

/** Forcing term for meso-micro scale coupled ABL simulations
 *  \ingroup icns_src we_abl
 *
 *  \sa ABL
 */

class ABLMesoForcingMom
    : public ABLMesoscaleForcing
    , public MomentumSource::Register<ABLMesoForcingMom>
{
public:
    static std::string identifier() { return "ABLMesoForcingMom"; }

    explicit ABLMesoForcingMom(const CFDSim& sim);

    ~ABLMesoForcingMom() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void mean_velocity_init(const ABLMesoscaleInput& ncfile);

    void mean_velocity_init(
        const VelPlaneAveragingFine& vavg, const ABLMesoscaleInput& ncfile);

    void
    mean_velocity_heights(std::unique_ptr<ABLMesoscaleInput> const& ncfile);

    void mean_velocity_heights(
        const VelPlaneAveragingFine& vavg,
        std::unique_ptr<ABLMesoscaleInput> const& ncfile);

    amrex::Vector<amrex::Real>& mom_u_error() { return m_err_U; }
    amrex::Vector<amrex::Real>& mom_v_error() { return m_err_V; }

private:
    // these are inputs
    amrex::Gpu::DeviceVector<amrex::Real> m_meso_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_meso_u_vals;
    amrex::Gpu::DeviceVector<amrex::Real> m_meso_v_vals;

    // these are the instantaneous planar averages
    amrex::Gpu::DeviceVector<amrex::Real> m_vavg_ht;

    // these specify the source term
    amrex::Gpu::DeviceVector<amrex::Real> m_error_meso_avg_U;
    amrex::Gpu::DeviceVector<amrex::Real> m_error_meso_avg_V;

    amrex::Vector<amrex::Real> m_err_U;
    amrex::Vector<amrex::Real> m_err_V;

    // these should be m_ind_polyOrder+1 (currently, cubic polynomials are
    // assumed)
    amrex::Array<amrex::Real, 4> m_poly_coeff_U;
    amrex::Array<amrex::Real, 4> m_poly_coeff_V;
};

} // namespace amr_wind::pde::icns

#endif // ABLMESOFORCINGMOM_H
