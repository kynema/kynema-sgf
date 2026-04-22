#ifndef HELICS_H
#define HELICS_H

#include "AMReX_AmrCore.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

#ifdef AMR_WIND_USE_HELICS
#include "helics/cpp98/CombinationFederate.hpp"
#include "helics/cpp98/helics.hpp"
#include "helics/cpp98/Federate.hpp"
#include "AMReX_REAL.H"

using namespace amrex::literals;
#endif

namespace amr_wind {

class HelicsStorage
{
public:
    explicit HelicsStorage(CFDSim& sim);
    ~HelicsStorage();

    void pre_advance_work();
    [[nodiscard]] bool is_activated() const { return m_helics_activated; }

#ifdef AMR_WIND_USE_HELICS
    HelicsTime m_currenttime = 0.0_rt;
    std::string m_fedinitstring = "--federates=1";
#endif

    int m_num_turbines{1};

    amrex::Real m_inflow_wind_speed_to_amrwind{8.0_rt};
    amrex::Real m_inflow_wind_direction_to_amrwind{270.0_rt};

    amrex::Vector<amrex::Real> m_turbine_power_to_controller;
    amrex::Vector<amrex::Real> m_turbine_wind_direction_to_controller;
    amrex::Vector<amrex::Real> m_turbine_yaw_to_amrwind;

private:
    CFDSim& m_sim;

    void send_messages_to_controller();
    void recv_messages_from_controller();

#ifdef AMR_WIND_USE_HELICS
    std::unique_ptr<helicscpp::FederateInfo> m_fi;
    std::unique_ptr<helicscpp::CombinationFederate> m_vfed;
#endif

#ifdef AMR_WIND_USE_HELICS
    bool m_helics_activated{true};
#else
    bool m_helics_activated{false};
#endif
};

} // namespace amr_wind

#endif /* HELICS_H */
