#ifndef ACTUATOR_H
#define ACTUATOR_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

/** \defgroup actuator Actuator source terms
 *  Actuator source modeling for wind turbines
 *  \ingroup wind
 */

namespace amr_wind {

class Field;

namespace actuator {

class ActuatorModel;
class ActuatorContainer;

/** Actuator line/disk modeling for wind turbines.
 *
 *  \ingroup actuator
 *
 *  This class provides an interface to model wind turbines in a wind farm as
 *  actuator line/disk source terms.
 *
 *  \sa ActuatorModel, ActuatorContainer
 */
class Actuator : public Physics::Register<Actuator>
{
public:
    static std::string identifier() { return "Actuator"; }

    explicit Actuator(CFDSim& /*sim*/);

    ~Actuator() override;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

    [[nodiscard]] ActuatorModel& get_act(int index) const
    {
        return *m_actuators.at(index);
    }

    [[nodiscard]] ActuatorModel&
    get_act_bylabel(const std::string& actlabel) const;

    template <typename T>
    T* get_actuator(std::string& key) const;

protected:
    //! Total number of actuator components (e.g., turbines) in the flow field
    [[nodiscard]] int num_actuators() const
    {
        return static_cast<int>(m_actuators.size());
    }

    virtual void prepare_outputs();

private:
    void setup_container();

    void update_positions();

    void update_velocities();

    void compute_forces();

    void compute_source_term();

    void communicate_turbine_io();

    CFDSim& m_sim;

    Field& m_act_source;

    std::vector<std::unique_ptr<ActuatorModel>> m_actuators;

    std::unique_ptr<ActuatorContainer> m_container;
};

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATOR_H */
