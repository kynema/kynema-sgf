#include "amr-wind/equation_systems/icns/source_terms/CoriolisForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Coriolis term for planetary rotation
 *
 *  Parameters are read from the `abl` namespace in the input file. The
 *  following parameters are available:
 *
 *  - `latitude`:  The latitude (in degrees) where the Coriolis term is
 * computed. This argument is mandatory.
 *
 *  - `east_vector`, `north_vector` (Optional) Vectors specifying the east and
 * north directions with respect to computational domain
 *
 *  - `rotational_time_period` Time period for planetary rotation (default:
 * 86400 seconds)
 */
CoriolisForcing::CoriolisForcing(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
{
    static_assert(AMREX_SPACEDIM == 3, "ABL implementation requires 3D domain");
    amrex::ParmParse pp("CoriolisForcing");

    // Latitude is mandatory, everything else is optional
    // Latitude is read in degrees
    pp.get("latitude", m_latitude);
    m_latitude = utils::radians(m_latitude);
    m_sinphi = std::sin(m_latitude);
    m_cosphi = std::cos(m_latitude);

    // Read the rotational time period (in seconds) -- This is 23hrs and 56
    // minutes and 4.091_rt seconds
    amrex::Real rot_time_period = 86164.091_rt;
    pp.query("rotational_time_period", rot_time_period);

    m_omega = utils::two_pi() / rot_time_period;
    m_coriolis_factor = 2.0_rt * m_omega;

    pp.queryarr("east_vector", m_east, 0, AMREX_SPACEDIM);
    pp.queryarr("north_vector", m_north, 0, AMREX_SPACEDIM);
    utils::vec_normalize(m_east.data());
    utils::vec_normalize(m_north.data());
    utils::cross_prod(m_east.data(), m_north.data(), m_up.data());

    pp.query("is_horizontal", m_is_horizontal);
}

CoriolisForcing::~CoriolisForcing() = default;

void CoriolisForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> east{
        m_east[0], m_east[1], m_east[2]};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> north{
        m_north[0], m_north[1], m_north[2]};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> up{m_up[0], m_up[1], m_up[2]};

    const auto sinphi = m_sinphi;
    const auto cosphi = m_cosphi;
    const auto corfac = m_coriolis_factor;
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);

    amrex::Real fac = (m_is_horizontal) ? 0. : 1.;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real ue = (east[0] * vel(i, j, k, 0)) +
                               (east[1] * vel(i, j, k, 1)) +
                               (east[2] * vel(i, j, k, 2));
        const amrex::Real un = (north[0] * vel(i, j, k, 0)) +
                               (north[1] * vel(i, j, k, 1)) +
                               (north[2] * vel(i, j, k, 2));
        const amrex::Real uu = (up[0] * vel(i, j, k, 0)) +
                               (up[1] * vel(i, j, k, 1)) +
                               (up[2] * vel(i, j, k, 2));

        const amrex::Real ae = +corfac * (un * sinphi - fac * uu * cosphi);
        const amrex::Real an = -corfac * ue * sinphi;
        const amrex::Real au = +fac * corfac * ue * cosphi;

        const amrex::Real ax = (ae * east[0]) + (an * north[0]) + (au * up[0]);
        const amrex::Real ay = (ae * east[1]) + (an * north[1]) + (au * up[1]);
        const amrex::Real az = (ae * east[2]) + (an * north[2]) + (au * up[2]);

        src_term(i, j, k, 0) += ax;
        src_term(i, j, k, 1) += ay;
        src_term(i, j, k, 2) += az;
    });
}

} // namespace amr_wind::pde::icns
