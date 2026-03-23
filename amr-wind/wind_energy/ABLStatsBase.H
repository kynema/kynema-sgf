#ifndef ABLSTATSBASE_H
#define ABLSTATSBASE_H

#include <cstdint>
#include "amr-wind/core/Factory.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLWallFunction.H"

namespace amr_wind {
namespace pde::icns {
class ABLForcing;
class ABLMesoForcingMom;
} // namespace pde::icns
namespace pde::temperature {
class ABLMesoForcingTemp;
}

enum class ABLStatsMode : std::uint8_t {
    computed = 0, ///! Computed using planar averages
    prescribed    ///! Prescribed as user inputs
};

class ABLStatsBase
    : public Factory<ABLStatsBase, CFDSim&, ABLWallFunction&, int>
{
public:
    static std::string base_identifier() { return "ABLStatsBase"; }

    ABLStatsBase() = default;

    ~ABLStatsBase() override = default;

    //! Flag indicating ABL simulation mode
    [[nodiscard]] virtual ABLStatsMode abl_mode() const = 0;

    //! Interpolating object for vertical velocity profile
    [[nodiscard]] virtual const VelPlaneAveraging&
    vel_profile_coarse() const = 0;
    [[nodiscard]] virtual const VelPlaneAveragingFine& vel_profile() const = 0;

    //! Interpolating object for vertical temperature profile
    [[nodiscard]] virtual const FieldPlaneAveraging& theta_profile() const = 0;
    [[nodiscard]] virtual const FieldPlaneAveragingFine&
    theta_profile_fine() const = 0;

    //! Perform initialization actions after the mesh has been created
    virtual void post_init_actions() = 0;

    //! Perform actions at the beginning of a timestep
    virtual void pre_advance_work() = 0;

    //! Perform actions at the end of a timestep
    virtual void post_advance_work() = 0;

    virtual void
    register_forcing_term(pde::icns::ABLForcing* forcing) const = 0;

    virtual void
    register_meso_mom_forcing(pde::icns::ABLMesoForcingMom* forcing) const = 0;

    virtual void register_meso_temp_forcing(
        pde::temperature::ABLMesoForcingTemp* forcing) const = 0;
};

} // namespace amr_wind

#endif /* ABLSTATSBASE_H */
