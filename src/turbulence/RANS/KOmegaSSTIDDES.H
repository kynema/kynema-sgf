#ifndef KOMEGASSTIDDES_H
#define KOMEGASSTIDDES_H

#include <string>
#include "amr-wind/turbulence/RANS/KOmegaSST.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::turbulence {

/** K-Omega-SST-IDDES Hybrid RANS/LES turbulence model
 *
 *
 * Turbulence model described in
 *
 * Gritskevich, M. S., Garbaruk, A.V., Schutze, J and Menter, F.R., "Development
 * of DDES and IDDES formulations for the k-omega Shear Stress Transport Model",
 * Flow Turbulence and Combustion, 88, 2012, pp. 431-449, DOI
 * 10.1007_rt/s10494-011-9378-4.
 *
 * Implementation contains the simplified version where f_e = 0 as in eq. 12 in
 * the paper
 *
 */
template <typename Transport>
class KOmegaSSTIDDES : public KOmegaSST<Transport>
{
public:
    // cppcheck-suppress duplInheritedMember
    static std::string identifier()
    {
        return "KOmegaSSTIDDES-" + Transport::identifier();
    }

    explicit KOmegaSSTIDDES(CFDSim& sim);

    ~KOmegaSSTIDDES() override;

    [[nodiscard]] std::string model_name() const override
    {
        return "KOmegaSSTIDDES";
    }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType diff_type) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

protected:
    Field& m_rans_ind;

    //! Turbulence constants
    amrex::Real m_Cdes1{0.78_rt};
    amrex::Real m_Cdes2{0.61_rt};
    amrex::Real m_Cdt1{20.0_rt};
    amrex::Real m_Cdt2{3.0_rt};
    amrex::Real m_Cl{5.0_rt};
    amrex::Real m_Ct{1.87_rt};
    amrex::Real m_Cw{0.15_rt};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_sdr_prod_clip_factor{10.0_rt};
};

} // namespace amr_wind::turbulence

#endif /* KOMEGASSTIDDES_H */
