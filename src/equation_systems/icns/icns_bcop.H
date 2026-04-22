#ifndef ICNS_BCOP_H
#define ICNS_BCOP_H

#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/BCOps.H"
#include "amr-wind/boundary_conditions/velocity_bcs.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/** BC operations for incompressible Navier-Stokes
 *  \ingroup icns
 */
template <>
struct BCOp<ICNS>
{
    BCOp(PDEFields& fields, const SimTime& time)
        : m_fields(fields), m_time(time)
    {}

    /** Initialize the BCs for the PDE variable as well as its source term
     */
    void init_bcs()
    {
        BCVelocity bc(m_fields.field);
        bc();
        const auto udfs = bc.get_dirichlet_udfs();
        vel_bc::register_velocity_dirichlet(
            m_fields.field, m_fields.repo.mesh(), m_time, udfs);

        auto& density = m_fields.repo.get_field("density");
        const amrex::Real density_default = 1.0_rt;
        BCScalar bc_den(density);
        bc_den(density_default);

        auto& pressure = m_fields.repo.get_field("p");
        BCPressure bc_pres(pressure);
        bc_pres();

        // Used for fillpatch operation on the source term
        BCSrcTerm bc_src(m_fields.src_term);
        bc_src();

        // Used for fillpatch operation during regrid
        auto& gradp = m_fields.repo.get_field("gp");
        gradp.set_default_fillpatch_bc(m_time);
    }

    /** Apply boundary conditions before a linear solve or compute_diff_term
     */
    void apply_bcs(const FieldState rho_state)
    {
        amrex::IntVect ng_diff(1);
        auto& field = m_fields.field;
        if ((rho_state != FieldState::Old && rho_state != FieldState::New)) {
            amrex::Abort(
                "icns_bcop.H apply_bcs(): a state other than New or Old was "
                "used. The time used for fillphysbc must be evaluated before "
                "using a different state for this routine.\n");
        }
        if (rho_state == FieldState::New) {
            field.fillphysbc(m_time.new_time(), ng_diff);
        }
        field.apply_bc_funcs(rho_state);
    }

    PDEFields& m_fields;
    const SimTime& m_time;
};

} // namespace amr_wind::pde

#endif /* ICNS_BCOP_H */
