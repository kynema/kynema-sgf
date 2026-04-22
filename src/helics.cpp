#include "amr-wind/helics.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"

#include <sstream>
#include <iostream>
#include <string>
#include <list>
#include "AMReX_REAL.H"

using namespace amrex::literals;

#ifdef AMR_WIND_USE_HELICS
using namespace helicscpp;
#endif

void tokenize(
    const std::string& s,
    const std::string& del,
    std::list<amrex::Real>& return_list)
{

    int start = 0;
    int end = static_cast<int>(s.find(del));

    return_list.push_front(
        static_cast<amrex::Real>(
            atof(s.substr(start + 1, end - start).c_str())));
    while (end > 0) {

        start = end + static_cast<int>(del.size());
        end = static_cast<int>(s.find(del, start));

        if (end == -1) {
            end = -2;
        }
        return_list.push_front(
            static_cast<amrex::Real>(
                atof(s.substr(start, end - start).c_str())));
    }
}

namespace amr_wind {

HelicsStorage::HelicsStorage(CFDSim& sim) : m_sim(sim)
{

#ifdef AMR_WIND_USE_HELICS
    amrex::ParmParse phelics("helics");
    phelics.query("activated", m_helics_activated);

    std::string broker_address = "127.0.0.1";
    int broker_port = 23405;

    phelics.query("broker_address", broker_address);
    phelics.query("broker_port", broker_port);

    m_fedinitstring.append(" --broker_address=" + broker_address);
    m_fedinitstring.append(" --brokerport=" + std::to_string(broker_port));
    amrex::Print() << "helics fedinitstring: " << m_fedinitstring << std::endl;

#endif

    if (!m_helics_activated) {
        return;
    }

#ifdef AMR_WIND_USE_HELICS
    if (amrex::ParallelDescriptor::IOProcessor()) {

        m_fi = std::make_unique<helicscpp::FederateInfo>("zmq");
        m_fi->setCoreInit(m_fedinitstring);

        m_vfed = std::make_unique<helicscpp::CombinationFederate>(
            "Test receiver Federate", *m_fi);

        m_vfed->setProperty(HELICS_PROPERTY_TIME_DELTA, 1.0_rt);
        //  Subscribe to PI SENDER's publication
        m_vfed->registerSubscription("control", "string");
        //  Register the publication
        m_vfed->registerGlobalPublication("status", "string");
        //  Enter initialization state
        m_vfed->enterInitializingMode(); // can throw
                                         // helicscpp::InvalidStateTransition
                                         // exception
        m_vfed->enterExecutingMode();
        amrex::Print() << "PI RECEIVER: Value federate created";
        amrex::Print() << "Creation complete!! ";
    }

    // parse input file and count how many actuators exist
    amrex::Vector<std::string> actuators;
    amrex::ParmParse pp("Actuator");
    pp.queryarr("labels", actuators);

    if (actuators.size() > 0) {
        m_num_turbines = actuators.size();
    }

    m_turbine_power_to_controller.resize(m_num_turbines, 0.0_rt);
    m_turbine_wind_direction_to_controller.resize(m_num_turbines, 0.0_rt);
    m_turbine_yaw_to_amrwind.resize(m_num_turbines, 270.0_rt);

#endif
}

void HelicsStorage::pre_advance_work()
{
#ifdef AMR_WIND_USE_HELICS
    if (m_helics_activated) {
        send_messages_to_controller();
        recv_messages_from_controller();
    }
#endif
}

void HelicsStorage::send_messages_to_controller()
{
#ifdef AMR_WIND_USE_HELICS
    if (amrex::ParallelDescriptor::IOProcessor()) {
        // put helics send stuff here
    }
#endif
}

void HelicsStorage::recv_messages_from_controller()
{

    amrex::Print() << "recv message from controller at time: "
                   << m_sim.time().current_time() << '\n';

#ifdef AMR_WIND_USE_HELICS
    // receive wind direction and speed from controller (1 + 1)
    // receive turbine yaw directions (num_turbines)
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::stringstream charFromControlCenter;
        m_currenttime = m_vfed->requestNextStep();
        int subCount = m_vfed->getInputCount();
        helicscpp::Input sub;

        for (int i = 0; i < subCount; i++) {
            int yy = i;
            sub = m_vfed->getSubscription(yy);
            charFromControlCenter << sub.getString().c_str();
            amrex::Print() << "\n advancing time "
                           << charFromControlCenter.str();
        }

        if (m_currenttime > 1) {
            // Igonre timestep 0 since message pipe has junk.
            // // Unpack the values from the control center using a string
            // stream

            std::list<amrex::Real> return_list;
            const std::string comma = ",";
            tokenize(charFromControlCenter.str(), comma, return_list);

            // skip the time stamp and source strings from return_list
            return_list.pop_front();
            return_list.pop_front();
            for (int i = m_turbine_yaw_to_amrwind.size() - 1; i >= 0; --i) {
                m_turbine_yaw_to_amrwind[i] = return_list.front();
                return_list.pop_front();
            }

            m_inflow_wind_direction_to_amrwind = return_list.front();
            return_list.pop_front();
            m_inflow_wind_speed_to_amrwind = return_list.front();
            return_list.pop_front();
            auto time = return_list.front();
            return_list.pop_front();

            amrex::Print() << "\n speed: " << m_inflow_wind_speed_to_amrwind
                           << "  direction: "
                           << m_inflow_wind_direction_to_amrwind << " Time "
                           << time << std::endl;
            for (int i = 0; i < m_turbine_yaw_to_amrwind.size(); ++i) {
                amrex::Print()
                    << "T" << i << " yaw: " << m_turbine_yaw_to_amrwind[i]
                    << ' ';
            }
            amrex::Print() << std::endl;
        }

        std::stringstream ssToControlCenter;

        helicscpp::Publication pub;

        int pubCount = m_vfed->getPublicationCount();

        amrex::Real wind_speed = 123.0_rt;
        amrex::Real wind_direction = 456.0_rt;

        auto& phy_mgr = m_sim.physics_manager();
        if (phy_mgr.contains("ABL")) {
            auto& abl = phy_mgr.get<amr_wind::ABL>();
            const amrex::Real height = 90.0_rt;
            wind_speed = abl.abl_statistics()
                             .vel_profile()
                             .line_hvelmag_average_interpolated(height);
            amrex::Real velx =
                abl.abl_statistics().vel_profile().line_average_interpolated(
                    height, 0);
            amrex::Real vely =
                abl.abl_statistics().vel_profile().line_average_interpolated(
                    height, 1);
            const amrex::Real turbine_angle = std::atan2(vely, velx);
            wind_direction =
                -amr_wind::utils::degrees(turbine_angle) + 270.0_rt;
        }

        amrex::Print() << "pub count: " << pubCount << std::endl;

        for (int i = 0; i < pubCount; i++) {
            pub = m_vfed->getPublication(i);

            ssToControlCenter << "[" << m_sim.time().current_time() << ", "
                              << wind_speed << " , " << wind_direction;

            for (int yy = 0; yy < m_turbine_power_to_controller.size(); yy++) {
                ssToControlCenter << ",";
                ssToControlCenter << m_turbine_power_to_controller[yy];
            }

            for (int yy = 0; yy < m_turbine_wind_direction_to_controller.size();
                 yy++) {
                ssToControlCenter << ",";
                ssToControlCenter << m_turbine_wind_direction_to_controller[yy];
            }

            ssToControlCenter << "]";

            std::string strToControlCenter;
            strToControlCenter = ssToControlCenter.str();
            pub.publish(strToControlCenter.c_str());

            amrex::Print() << "\n should send m_turbine_power_to_controller "
                           << strToControlCenter.c_str() << std::endl;
        }
    }

    // broadcast wind turbine yaw angles to all procs
    // TODO: some day only need to send/recv to specific turbines
    amrex::ParallelDescriptor::Bcast(
        m_turbine_yaw_to_amrwind.data(), m_turbine_yaw_to_amrwind.size(),
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());

    // broadcast wind speed and direction to all procs
    amrex::ParallelDescriptor::Bcast(
        &m_inflow_wind_speed_to_amrwind, 1,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());

    amrex::ParallelDescriptor::Bcast(
        &m_inflow_wind_direction_to_amrwind, 1,
        amrex::ParallelDescriptor::IOProcessorNumber(),
        amrex::ParallelDescriptor::Communicator());
#endif
}

HelicsStorage::~HelicsStorage() = default;

} // namespace amr_wind
