#ifndef ABLANELASTIC_H
#define ABLANELASTIC_H
#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/utilities/constants.H"
#include "amr-wind/utilities/MultiLevelVector.H"
#include "amr-wind/eos_models/EOSModel.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/** Interface for ABL anelastic
 *  \ingroup we_abl
 *
 *  This class performs the necessary operations for anelastic ABL
 *
 */
class ABLAnelastic
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    explicit ABLAnelastic(CFDSim& /*sim*/);

    void post_init_actions();
    void post_regrid_actions();
    void initialize_data();
    [[nodiscard]] bool is_anelastic() const { return m_is_anelastic; }

private:
    //! Function to calculate the hydrostatic density and pressure
    void initialize_isentropic_hse();

    const CFDSim& m_sim;

    bool m_is_anelastic{false};

    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    amrex::Real m_reference_density_constant{1.0_rt};

    amrex::Real m_bottom_reference_pressure{1.0e5_rt};

    const int m_axis{2};

    MultiLevelVector m_density{FieldLoc::CELL};
    MultiLevelVector m_pressure{FieldLoc::NODE};
    MultiLevelVector m_theta{FieldLoc::CELL};
};
} // namespace amr_wind
#endif
