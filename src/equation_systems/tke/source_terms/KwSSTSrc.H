#ifndef KWSSTSRC_H
#define KWSSTSRC_H

#include "src/equation_systems/tke/TKESource.H"

namespace kynema_sgf::pde::tke {

/** Source terms for TKE equation in K-w-SST turbulence model
 *   \ingroup turb_model
 *
 *  Works for all TKE equations with production and dissipation terms.
 *
 *  For FieldState::Old (the Godunov forcing) the source is the full explicit
 *  balance. For the other states it also adds back the part of the
 *  dissipation that the turbulence model linearized onto the diagonal of the
 *  linear system (`tke_lhs_src_term`), so that a converged solution satisfies
 *  the full balance independent of the diffusion type and time step.
 *
 *  The model terms are density weighted; the source is returned per unit mass
 *  because the TKE equation multiplies it by density.
 */
class KwSSTSrc : public TKESource::Register<KwSSTSrc>
{
public:
    static std::string identifier() { return "KwSSTSrc"; }

    explicit KwSSTSrc(const CFDSim& sim);

    ~KwSSTSrc() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    Field& m_shear_prod;
    Field& m_diss;
    Field& m_buoy_term;
    Field& m_tke;
    Field& m_lhs;
    Field& m_density;
};

} // namespace kynema_sgf::pde::tke

#endif /* KWSSTSRC_H */
