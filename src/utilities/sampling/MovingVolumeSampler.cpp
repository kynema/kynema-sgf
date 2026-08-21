#include "src/utilities/sampling/MovingVolumeSampler.H"

#include "src/CFDSim.H"
#include "src/utilities/index_operations.H"
#include "src/utilities/sampling/SamplingUtils.H"
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
    if (lhs.size() != rhs.size()) {
        return true;
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (amrex::Math::abs(lhs[d] - rhs[d]) > constants::TIGHT_TOL) {
            return true;
        }
    }
    return false;
}

} // namespace

MovingVolumeSampler::MovingVolumeSampler(const CFDSim& sim) : VolumeSampler(sim)
{}

void MovingVolumeSampler::initialize(const std::string& key)
{
    VolumeSampler::initialize(key);
    m_reference.origin = m_lo;
    m_reference.axis1 = {m_hi[0] - m_lo[0], 0.0, 0.0};
    m_reference.axis2 = {0.0, m_hi[1] - m_lo[1], 0.0};
    m_reference.axis3 = {0.0, 0.0, m_hi[2] - m_lo[2]};

    amrex::ParmParse pp(key);
    pp.get("actuator_label", m_actuator_label);
    if (!m_sim.physics_manager().contains(actuator::Actuator::identifier())) {
        amrex::Abort(key + " requires Actuator in incflo.physics");
    }

    const auto& actuators = m_sim.physics_manager().get<actuator::Actuator>();
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

void MovingVolumeSampler::check_bounds()
{
    if (m_defer_bounds_check) {
        return;
    }

    const auto* prob_lo = m_sim.mesh().Geom(0).ProbLo();
    const auto* prob_hi = m_sim.mesh().Geom(0).ProbHi();
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const auto point =
                        m_geometry.origin[d] + (i * m_geometry.axis1[d]) +
                        (j * m_geometry.axis2[d]) + (k * m_geometry.axis3[d]);
                    if ((point < prob_lo[d]) || (point >= prob_hi[d])) {
                        amrex::Abort(
                            "MovingVolumeSampler: Point out of domain. "
                            "Redefine the volume so it remains completely "
                            "inside the domain.");
                    }
                }
            }
        }
    }
}

bool MovingVolumeSampler::update_geometry(const amrex::Real time)
{
    AMREX_ALWAYS_ASSERT(m_actuator != nullptr);
    const auto frame = m_actuator->reference_frame(time);
    if (!frame.has_value()) {
        amrex::Abort(
            "MovingVolumeSampler: actuator '" + m_actuator_label +
            "' does not provide a rigid reference frame");
        return false;
    }

    const auto& rigid_frame = *frame;
    auto origin = values(rigid_frame.apply_point(vector(m_reference.origin)));
    auto axis1 = values(rigid_frame.apply_vector(vector(m_reference.axis1)));
    auto axis2 = values(rigid_frame.apply_vector(vector(m_reference.axis2)));
    auto axis3 = values(rigid_frame.apply_vector(vector(m_reference.axis3)));
    const bool geometry_changed = changed(origin, m_geometry.origin) ||
                                  changed(axis1, m_geometry.axis1) ||
                                  changed(axis2, m_geometry.axis2) ||
                                  changed(axis3, m_geometry.axis3);

    m_geometry = {
        std::move(origin), std::move(axis1), std::move(axis2),
        std::move(axis3)};
    check_bounds();
    return geometry_changed;
}

bool MovingVolumeSampler::update_sampling_locations()
{
    return update_geometry(m_sim.time().new_time());
}

void MovingVolumeSampler::sampling_locations(SampleLocType& sample_locs) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());
    sampling_locations(sample_locs, m_sim.mesh().Geom(0).Domain());
    AMREX_ALWAYS_ASSERT(sample_locs.locations().size() == num_points());
}

void MovingVolumeSampler::sampling_locations(
    SampleLocType& sample_locs, const amrex::Box& box) const
{
    AMREX_ALWAYS_ASSERT(sample_locs.locations().empty());

    const auto& dxinv = m_sim.mesh().Geom(0).InvCellSizeArray();
    const auto& plo = m_sim.mesh().Geom(0).ProbLoArray();
    const auto& fine_geom = m_sim.mesh().Geom(m_sim.mesh().finestLevel());
    int idx = 0;
    for (int k = 0; k < m_npts_dir[2]; ++k) {
        for (int j = 0; j < m_npts_dir[1]; ++j) {
            for (int i = 0; i < m_npts_dir[0]; ++i) {
                amrex::RealVect loc;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    loc[d] = m_geometry.origin[d] +
                             (m_geometry.axis1[d] * i / m_npts_dir[0]) +
                             (m_geometry.axis2[d] * j / m_npts_dir[1]) +
                             (m_geometry.axis3[d] * k / m_npts_dir[2]);
                }

                if (m_snap_to_cell_center) {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        loc[d] = sampling_utils::snap_to_nearest_cell_center(
                            fine_geom, d, loc[d]);
                    }
                }

                if (utils::contains(box, loc, plo, dxinv)) {
                    sample_locs.push_back(loc, idx);
                }
                ++idx;
            }
        }
    }
}

#ifdef KYNEMA_SGF_USE_NETCDF
void MovingVolumeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& grp) const
{
    VolumeSampler::define_netcdf_metadata(grp);
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("actuator_label", m_actuator_label);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void MovingVolumeSampler::output_netcdf_data(
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
void MovingVolumeSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}

void MovingVolumeSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}
#endif

void MovingVolumeSampler::populate_info_file(std::ostream& fh) const
{
    fh << "   actuator_label: " << m_actuator_label << "\n";
}

} // namespace sampling

template struct ::kynema_sgf::sampling::SamplerBase::Register<
    ::kynema_sgf::sampling::MovingVolumeSampler>;

} // namespace kynema_sgf
