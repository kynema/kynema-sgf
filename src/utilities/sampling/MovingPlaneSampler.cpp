#include "src/utilities/sampling/MovingPlaneSampler.H"

#include "src/CFDSim.H"
#include "src/wind_energy/actuator/Actuator.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

#include "AMReX_ParmParse.H"

namespace kynema_sgf {
namespace sampling {
namespace {

amrex::Vector<amrex::Real> values(const vs::Vector& vector)
{
    return {vector.x(), vector.y(), vector.z()};
}

vs::Vector vector(const amrex::Vector<amrex::Real>& values)
{
    return {values[0], values[1], values[2]};
}

bool changed(
    const amrex::Vector<amrex::Real>& lhs,
    const amrex::Vector<amrex::Real>& rhs)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (amrex::Math::abs(lhs[d] - rhs[d]) > constants::TIGHT_TOL) {
            return true;
        }
    }
    return false;
}

} // namespace

MovingPlaneSampler::MovingPlaneSampler(const CFDSim& sim) : PlaneSampler(sim)
{}

void MovingPlaneSampler::initialize(const std::string& key)
{
    PlaneSampler::initialize(key);
    m_reference = {m_origin, m_axis1, m_axis2, m_offset_vector};

    amrex::ParmParse pp(key);
    pp.get("actuator_label", m_actuator_label);
    if (!m_sim.physics_manager().contains(actuator::Actuator::identifier())) {
        amrex::Abort(key + " requires Actuator in incflo.physics");
    }

    const auto& actuators =
        m_sim.physics_manager().get<actuator::Actuator>();
    const auto& model = actuators.get_act_bylabel(m_actuator_label);
    if (model.label() != m_actuator_label) {
        amrex::Abort(
            key + ": cannot find actuator with label '" + m_actuator_label +
            "'");
    }
    m_actuator = &model;
    m_defer_bounds_check = false;
    update_geometry(m_sim.time().current_time());
}

void MovingPlaneSampler::check_bounds()
{
    if (!m_defer_bounds_check) {
        PlaneSampler::check_bounds();
    }
}

bool MovingPlaneSampler::update_geometry(const amrex::Real time)
{
    AMREX_ALWAYS_ASSERT(m_actuator != nullptr);
    const auto frame = m_actuator->reference_frame(time);
    if (!frame.has_value()) {
        amrex::Abort(
            "MovingPlaneSampler: actuator '" + m_actuator_label +
            "' does not provide a rigid reference frame");
    }

    auto origin = values(frame->apply_point(vector(m_reference.origin)));
    auto axis1 = values(frame->apply_vector(vector(m_reference.axis1)));
    auto axis2 = values(frame->apply_vector(vector(m_reference.axis2)));
    auto offset_vector =
        values(frame->apply_vector(vector(m_reference.offset_vector)));
    const bool geometry_changed =
        changed(origin, m_origin) || changed(axis1, m_axis1) ||
        changed(axis2, m_axis2) ||
        changed(offset_vector, m_offset_vector);

    m_origin = std::move(origin);
    m_axis1 = std::move(axis1);
    m_axis2 = std::move(axis2);
    m_offset_vector = std::move(offset_vector);
    PlaneSampler::check_bounds();
    return geometry_changed;
}

bool MovingPlaneSampler::update_sampling_locations()
{
    return update_geometry(m_sim.time().new_time());
}

#ifdef KYNEMA_SGF_USE_NETCDF
void MovingPlaneSampler::define_netcdf_metadata(
    const ncutils::NCGroup& grp) const
{
    PlaneSampler::define_netcdf_metadata(grp);
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("actuator_label", m_actuator_label);
    grp.def_var(
        "points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void MovingPlaneSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    SampleLocType sample_locs;
    sampling_locations(sample_locs);
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{
        1, static_cast<size_t>(num_points()), AMREX_SPACEDIM};
    grp.var("points").put(sample_locs.locations()[0].begin(), start, count);
}
#else
void MovingPlaneSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void MovingPlaneSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}
#endif

void MovingPlaneSampler::populate_info_file(std::ostream& fh) const
{
    PlaneSampler::populate_info_file(fh);
    fh << "   actuator_label: " << m_actuator_label << "\n";
}

} // namespace sampling

template struct ::kynema_sgf::sampling::SamplerBase::Register<
    ::kynema_sgf::sampling::MovingPlaneSampler>;

} // namespace kynema_sgf
