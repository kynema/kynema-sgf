#ifndef EOSMODEL_H
#define EOSMODEL_H

#include "AMReX_REAL.H"
#include "amr-wind/utilities/constants.H"

using namespace amrex::literals;

namespace amr_wind::eos {

/**
 *  \defgroup eos Equation of state models
 *
 *  AMR-Wind representation of equation of state models. These
 *  expression are adapted from those found in ERF
 *  https://github.com/erf-model/ERF.
 *
 */

/** Gamma law equation of state
 *  \ingroup eos
 */
struct GammaLaw
{
    using eos_type = GammaLaw;
    static std::string identifier() { return "GammaLaw"; }

    /**
     * Return pressure given density and potential temperature
     *
     * \param [in] rho    density
     * \param [in] theta  potential temperature
     * \param [in] qv     water vapor
     * \return            pressure
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real p_rth(
        const amrex::Real rho,
        const amrex::Real theta,
        const amrex::Real qv = 0.0_rt) const
    {
        const auto rhotheta = rho * theta;
        return m_p0 * std::pow(
                          m_air_gas_constant * rhotheta *
                              (1.0_rt + (m_water_vapor_gas_constant /
                                         m_air_gas_constant) *
                                            qv) *
                              m_ip0,
                          m_gamma);
    }

    /**
     * Return dP/drho at constant theta
     *
     * \param [in] rho   density
     * \param [in] theta potential temperature
     * \param [in] qv    water vapor
     * \return           pressure
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real dp_constanttheta(
        const amrex::Real rho,
        const amrex::Real theta,
        const amrex::Real qv = 0.0_rt) const
    {
        return m_gamma * m_p0 *
               std::pow(
                   (m_air_gas_constant * theta *
                    (1.0_rt +
                     m_water_vapor_gas_constant / m_air_gas_constant * qv) *
                    m_ip0),
                   m_gamma) *
               std::pow(rho, m_gamma - 1.0_rt);
    }

    template <class... Args>
    AMREX_GPU_HOST_DEVICE explicit GammaLaw(const amrex::Real p0 = 1.0e-5_rt)
        : m_p0{p0}, m_ip0{1.0_rt / p0}
    {}

    const amrex::Real m_gamma{constants::HEAT_CAPACITY_RATIO};
    const amrex::Real m_p0{1.0e5_rt};
    const amrex::Real m_ip0{1.0e-5_rt};
    const amrex::Real m_air_gas_constant{
        constants::UNIVERSAL_GAS_CONSTANT / constants::MOLAR_MASS_AIR};
    const amrex::Real m_water_vapor_gas_constant{
        constants::UNIVERSAL_GAS_CONSTANT / constants::MOLAR_MASS_WATER_VAPOR};
};
} // namespace amr_wind::eos

#endif /* EOSMODEL_H */
