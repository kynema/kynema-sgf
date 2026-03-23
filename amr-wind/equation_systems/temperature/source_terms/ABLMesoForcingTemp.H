#ifndef ABLMESOFORCINGTEMP_H
#define ABLMESOFORCINGTEMP_H

#include <memory>
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/wind_energy/ABLMesoscaleForcing.H"
#include "amr-wind/wind_energy/ABLMesoscaleInput.H"
#include "AMReX_REAL.H"

namespace amr_wind::pde::temperature {

class ABLMesoForcingTemp
    : public ABLMesoscaleForcing
    , public TemperatureSource::Register<ABLMesoForcingTemp>
{
public:
    static std::string identifier() { return "ABLMesoForcingTemp"; }

    explicit ABLMesoForcingTemp(const CFDSim& sim);

    ~ABLMesoForcingTemp() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void mean_temperature_init(const ABLMesoscaleInput& ncfile);

    void mean_temperature_init(
        const FieldPlaneAveragingFine& tavg, const ABLMesoscaleInput& ncfile);

    amrex::Real
    mean_temperature_heights(std::unique_ptr<ABLMesoscaleInput> const& ncfile);

    amrex::Real mean_temperature_heights(
        const FieldPlaneAveragingFine& tavg,
        std::unique_ptr<ABLMesoscaleInput> const& ncfile);

    amrex::Vector<amrex::Real>& theta_error() { return m_err_Theta; }

private:
    // these are inputs
    amrex::Gpu::DeviceVector<amrex::Real> m_meso_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_meso_theta_vals;

    // this is the instantaneous planar average (at AMR-Wind levels)
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_ht;

    // this specifies the source term
    amrex::Gpu::DeviceVector<amrex::Real> m_error_meso_avg_theta;

    amrex::Vector<amrex::Real> m_err_Theta;

    // this should be m_ind_polyOrder+1 (currently, cubic polynomial is assumed)
    amrex::Array<amrex::Real, 4> m_poly_coeff_theta;
};

} // namespace amr_wind::pde::temperature

#endif
