#ifndef KSGSM84SRC_H
#define KSGSM84SRC_H

#include "amr-wind/equation_systems/tke/TKESource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::tke {

/** TKE source term based on Moeng 1984 paper
 *  \ingroup tke_src turb_model we_abl
 */
class KsgsM84Src : public TKESource::Register<KsgsM84Src>
{
public:
    static std::string identifier() { return "KsgsM84Src"; }

    explicit KsgsM84Src(const CFDSim& /*sim*/);

    ~KsgsM84Src() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    amrex::Real m_Ceps{0.93_rt};
    amrex::Real m_CepsGround{3.9_rt};

    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_tke;
};

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real calc_ceps_local(
    const amrex::Real Ceps, const amrex::Real tlscale, const amrex::Real ds)
{
    return (Ceps / 0.93_rt) * (0.19_rt + (0.74_rt * tlscale / ds));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real calc_dissip(
    const amrex::Real Ceps, const amrex::Real tke, const amrex::Real tlscale)
{
    return Ceps * std::sqrt(tke) * tke / tlscale;
}

} // namespace amr_wind::pde::tke

#endif /* KSGSM84SRC_H */
