#ifndef SHEARSTRESSSIMPLE_H
#define SHEARSTRESSSIMPLE_H

#include "amr-wind/boundary_conditions/wall_models/LogLaw.H"
#include "amr-wind/wind_energy/ShearStress.H"
#include "amr-wind/boundary_conditions/wall_models/MOSD.H"

namespace amr_wind {

struct SimpleShearSchumann
{
    explicit SimpleShearSchumann(const amr_wind::LogLaw& ll)
        : utau2(ll.utau_mean * ll.utau_mean), wspd_mean(ll.wspd_mean), m_ll(ll)
    {}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real get_shear(
        amrex::Real u,
        amrex::Real /* wspd */,
        amrex::Real /*u_dx*/,
        amrex::Real /*v_dx*/,
        amrex::Real /*x_c*/,
        amrex::Real /*unit_nor*/) const
    {
        return u / wspd_mean * utau2;
    };

    amrex::Real utau2;
    amrex::Real wspd_mean;
    const amr_wind::LogLaw m_ll;
};
struct SimpleShearLogLaw
{
    explicit SimpleShearLogLaw(const amr_wind::LogLaw& ll) : m_ll(ll) {}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real get_shear(
        amrex::Real u,
        amrex::Real wspd,
        amrex::Real /*u_dx*/,
        amrex::Real /*v_dx*/,
        amrex::Real /*x_c*/,
        amrex::Real /*unit_nor*/) const
    {
        amrex::Real utau = m_ll.get_utau(wspd);
        return utau * utau * u / wspd;
    };

    const amr_wind::LogLaw m_ll;
};

struct SimpleShearMOSD
{
    explicit SimpleShearMOSD(
        const amr_wind::LogLaw& ll, const amr_wind::MOSD& md)
        : m_ll(ll), m_md(md)
    {}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real get_shear(
        amrex::Real u,
        amrex::Real wspd,
        amrex::Real u_dx,
        amrex::Real v_dx,
        amrex::Real x_c,
        amrex::Real unit_nor) const
    {
        amrex::Real utau = m_ll.get_utau(wspd);
        amrex::Real tau_vis = utau * utau * u / wspd;

        amrex::Real tau_wave = m_md.get_dyn_tau(u_dx, v_dx, x_c, unit_nor);

        return tau_vis + tau_wave;
    };

    const amr_wind::LogLaw m_ll;
    const amr_wind::MOSD m_md;
};

} // namespace amr_wind

#endif /* SHEARSTRESSSIMPLE_H */
