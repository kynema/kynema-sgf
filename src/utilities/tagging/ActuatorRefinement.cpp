#include "src/utilities/tagging/ActuatorRefinement.H"

#include "src/CFDSim.H"
#include "src/utilities/index_operations.H"
#include "src/wind_energy/actuator/Actuator.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

#include <algorithm>
#include <cmath>

#include "AMReX.H"
#include "AMReX_ParmParse.H"

using namespace amrex::literals;

namespace kynema_sgf {
namespace {

struct RefinementCylinder
{
    vs::Vector center;
    vs::Vector normal;
    amrex::Real radius;
    amrex::Real forward;
    amrex::Real backward;
    amrex::RealBox bound_box;
};

amrex::RealBox cylinder_bound_box(
    const vs::Vector& center,
    const vs::Vector& normal,
    const amrex::Real radius,
    const amrex::Real forward,
    const amrex::Real backward)
{
    const auto start = center - backward * normal;
    const auto end = center + forward * normal;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> lo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> hi;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        const amrex::Real radial_extent =
            radius *
            std::sqrt(amrex::max(0.0_rt, 1.0_rt - normal[n] * normal[n]));
        lo[n] = amrex::min(start[n], end[n]) - radial_extent;
        hi[n] = amrex::max(start[n], end[n]) + radial_extent;
    }
    return {lo.data(), hi.data()};
}

void require_nonnegative(
    const std::string& key, const std::string& name, const amrex::Real value)
{
    if (value < 0.0_rt) {
        amrex::Abort(
            key + "." + name + " must be greater than or equal to zero");
    }
}

bool contains_any(
    amrex::ParmParse& pp, const std::initializer_list<const char*>& names)
{
    return std::any_of(names.begin(), names.end(), [&pp](const char* name) {
        return pp.contains(name);
    });
}

} // namespace

ActuatorRefinement::ActuatorRefinement(CFDSim& sim)
    : m_sim(sim), m_max_level(sim.mesh().maxLevel() - 1)
{}

void ActuatorRefinement::initialize(const std::string& key)
{
    m_key = key;
    auto& physics = m_sim.physics_manager();
    if (!physics.contains(actuator::Actuator::identifier())) {
        amrex::Abort(key + " requires Actuator in incflo.physics");
    }

    amrex::ParmParse pp(key);
    pp.getarr("actuator_labels", m_actuator_labels);
    if (m_actuator_labels.empty()) {
        amrex::Abort(key + ".actuator_labels must not be empty");
    }

    pp.query("min_level", m_min_level);
    pp.query("max_level", m_max_level);
    if (m_min_level < 0 || m_max_level < m_min_level ||
        m_max_level >= m_sim.mesh().maxLevel()) {
        amrex::Abort(
            key + " requires 0 <= min_level <= max_level < amr.max_level");
    }

    pp.query("radial_padding_epsilon", m_radial_padding_epsilon);
    pp.query("radial_padding", m_radial_padding);

    const bool has_epsilon = contains_any(
        pp, {"axial_padding_epsilon", "forward_padding_epsilon",
             "backward_padding_epsilon"});
    const bool has_diameter = contains_any(
        pp, {"axial_padding_diameter", "forward_padding_diameter",
             "backward_padding_diameter"});
    const bool has_absolute = contains_any(
        pp, {"axial_padding", "forward_padding", "backward_padding"});
    if (static_cast<int>(has_epsilon) + static_cast<int>(has_diameter) +
            static_cast<int>(has_absolute) >
        1) {
        amrex::Abort(
            key +
            " must specify axial padding in only one of epsilon, diameter, "
            "or absolute units");
    }

    if (has_diameter) {
        m_axial_padding_type = AxialPaddingType::Diameter;
        m_forward_axial_padding = 0.0_rt;
        m_backward_axial_padding = 0.0_rt;
        amrex::Real shared = 0.0_rt;
        if (pp.query("axial_padding_diameter", shared)) {
            m_forward_axial_padding = shared;
            m_backward_axial_padding = shared;
        }
        pp.query("forward_padding_diameter", m_forward_axial_padding);
        pp.query("backward_padding_diameter", m_backward_axial_padding);
    } else if (has_absolute) {
        m_axial_padding_type = AxialPaddingType::Absolute;
        m_forward_axial_padding = 0.0_rt;
        m_backward_axial_padding = 0.0_rt;
        amrex::Real shared = 0.0_rt;
        if (pp.query("axial_padding", shared)) {
            m_forward_axial_padding = shared;
            m_backward_axial_padding = shared;
        }
        pp.query("forward_padding", m_forward_axial_padding);
        pp.query("backward_padding", m_backward_axial_padding);
    } else {
        m_axial_padding_type = AxialPaddingType::Epsilon;
        amrex::Real shared = 3.0_rt;
        pp.query("axial_padding_epsilon", shared);
        m_forward_axial_padding = shared;
        m_backward_axial_padding = shared;
        pp.query("forward_padding_epsilon", m_forward_axial_padding);
        pp.query("backward_padding_epsilon", m_backward_axial_padding);
    }

    require_nonnegative(
        key, "radial_padding_epsilon", m_radial_padding_epsilon);
    require_nonnegative(key, "radial_padding", m_radial_padding);
    require_nonnegative(key, "forward axial padding", m_forward_axial_padding);
    require_nonnegative(
        key, "backward axial padding", m_backward_axial_padding);
}

void ActuatorRefinement::resolve_actuators()
{
    const auto& actuators = m_sim.physics_manager().get<actuator::Actuator>();
    m_actuators.reserve(m_actuator_labels.size());
    for (const auto& label : m_actuator_labels) {
        const auto& model = actuators.get_act_bylabel(label);
        if (model.label() != label) {
            amrex::Abort(
                m_key + ": cannot find actuator with label '" + label + "'");
        }
        const auto geometries = model.refinement_geometries(0.0_rt);
        if (geometries.empty()) {
            amrex::Abort(
                m_key + ": actuator '" + label +
                "' does not provide rotor refinement geometry");
        }
        m_actuators.push_back(&model);
    }
}

void ActuatorRefinement::operator()(
    const int level,
    amrex::TagBoxArray& tags,
    const amrex::Real time,
    const int /*ngrow*/)
{
    if (level < m_min_level || level > m_max_level) {
        return;
    }
    if (m_actuators.empty()) {
        resolve_actuators();
    }

    amrex::Vector<RefinementCylinder> cylinders;
    for (const auto* actuator : m_actuators) {
        for (const auto& geometry : actuator->refinement_geometries(time)) {
            const amrex::Real radius =
                geometry.rotor_radius +
                m_radial_padding_epsilon * geometry.epsilon_max +
                m_radial_padding;
            amrex::Real axial_scale = 1.0_rt;
            if (m_axial_padding_type == AxialPaddingType::Epsilon) {
                axial_scale = geometry.epsilon_max;
            } else if (m_axial_padding_type == AxialPaddingType::Diameter) {
                axial_scale = 2.0_rt * geometry.rotor_radius;
            }
            const amrex::Real forward = m_forward_axial_padding * axial_scale;
            const amrex::Real backward = m_backward_axial_padding * axial_scale;
            cylinders.push_back(
                {geometry.center, geometry.normal, radius, forward, backward,
                 cylinder_bound_box(
                     geometry.center, geometry.normal, radius, forward,
                     backward)});
        }
    }

    const auto& geom = m_sim.mesh().Geom(level);
    const auto& field_fab = (*m_sim.repo().fields()[0])(level);
    for (amrex::MFIter mfi(field_fab); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tag = tags.array(mfi);
        for (const auto& cylinder : cylinders) {
            const auto tagged_box =
                bx & utils::realbox_to_box(cylinder.bound_box, geom);
            if (tagged_box.isEmpty()) {
                continue;
            }

            const auto center = cylinder.center;
            const auto normal = cylinder.normal;
            const amrex::Real radius_sq = cylinder.radius * cylinder.radius;
            const amrex::Real forward = cylinder.forward;
            const amrex::Real backward = cylinder.backward;
            const auto& problo = geom.ProbLoArray();
            const auto& dx = geom.CellSizeArray();

            amrex::ParallelFor(
                tagged_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    const vs::Vector point(
                        problo[0] + (i + 0.5_rt) * dx[0],
                        problo[1] + (j + 0.5_rt) * dx[1],
                        problo[2] + (k + 0.5_rt) * dx[2]);
                    const auto relative = point - center;
                    const amrex::Real axial = relative & normal;
                    const amrex::Real radial_sq = amrex::max(
                        0.0_rt, (relative & relative) - axial * axial);
                    if (axial <= forward && axial >= -backward &&
                        radial_sq <= radius_sq) {
                        tag(i, j, k) = amrex::TagBox::SET;
                    }
                });
        }
    }
}

} // namespace kynema_sgf
