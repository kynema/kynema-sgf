#ifndef ACTUATORREFINEMENT_H
#define ACTUATORREFINEMENT_H

#include "src/utilities/tagging/RefinementCriteria.H"
#include "src/wind_energy/actuator/actuator_types.H"

#include <string>

namespace kynema_sgf {

class CFDSim;

namespace actuator {
class ActuatorModel;
}

/** AMR refinement in moving cylinders surrounding actuator rotors. */
class ActuatorRefinement
    : public RefinementCriteria::Register<ActuatorRefinement>
{
public:
    static std::string identifier() { return "ActuatorRefinement"; }

    explicit ActuatorRefinement(CFDSim& sim);

    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    enum class AxialPaddingType { Epsilon, Diameter, Absolute };

    void resolve_actuators();

    CFDSim& m_sim;
    std::string m_key;
    amrex::Vector<std::string> m_actuator_labels;
    amrex::Vector<const actuator::ActuatorModel*> m_actuators;

    amrex::Real m_radial_padding_epsilon{3.0};
    amrex::Real m_radial_padding{0.0};
    amrex::Real m_forward_axial_padding{3.0};
    amrex::Real m_backward_axial_padding{3.0};
    AxialPaddingType m_axial_padding_type{AxialPaddingType::Epsilon};

    int m_min_level{0};
    int m_max_level{-1};
};

} // namespace kynema_sgf

#endif
