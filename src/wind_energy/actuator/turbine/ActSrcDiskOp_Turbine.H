#ifndef ACTSRCDISKOP_TURBINE_H
#define ACTSRCDISKOP_TURBINE_H

#include <numbers>
#include "AMReX_REAL.H"
#include "src/wind_energy/actuator/actuator_ops.H"
#include "src/wind_energy/actuator/actuator_utils.H"
#include "src/core/FieldRepo.H"
#include "src/wind_energy/actuator/turbine/turbine_types.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::ops {

template <typename ActTrait>
class ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<TurbineType, ActTrait>>>
{
private:
    typename ActTrait::DataType& m_data;
    Field& m_act_src;

    DeviceVecList m_pos;
    DeviceVecList m_force;
    DeviceVecList m_epsilon;
    DeviceTensorList m_orientation;

    DeviceVecComponent m_blades;
    DeviceVecComponent m_tower;
    DeviceVecComponent m_hub;

    void copy_to_device();

    //! Map a component view that points into host grid data onto the
    //! equivalent locations in the device copies of the grid data
    ComponentView to_device_view(const ComponentView& host_cv);

public:
    explicit ActSrcOp(typename ActTrait::DataType& data)
        : m_data(data)
        , m_act_src(m_data.sim().repo().get_field("actuator_src_term"))
    {}

    void initialize();

    void setup_op() { copy_to_device(); }

    void
    operator()(int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom);
};

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<TurbineType, ActTrait>>>::initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_force.resize(grid.force.size());
    m_epsilon.resize(grid.epsilon.size());
    m_orientation.resize(grid.orientation.size());

    const auto& meta = m_data.meta();
    m_blades.resize(meta.num_blades);
    m_tower.resize(1);
    m_hub.resize(1);
}

namespace detail {
template <typename T, typename HostVec, typename DevVec>
::kynema_sgf::utils::Slice<T> rebase_slice(
    const ::kynema_sgf::utils::Slice<T>& host_slice, const HostVec& host, DevVec& dev)
{
    if (host_slice.size() == 0) {
        return ::kynema_sgf::utils::Slice<T>{};
    }
    const auto offset = host_slice.data() - host.data();
    AMREX_ALWAYS_ASSERT(
        (offset >= 0) &&
        (static_cast<std::size_t>(offset) + host_slice.size() <=
         static_cast<std::size_t>(dev.size())));
    return ::kynema_sgf::utils::Slice<T>{dev.data() + offset, host_slice.size()};
}
} // namespace detail

template <typename ActTrait>
ComponentView ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<TurbineType, ActTrait>>>::
    to_device_view(const ComponentView& host_cv)
{
    const auto& grid = m_data.grid();
    ComponentView cv;
    cv.pos = detail::rebase_slice(host_cv.pos, grid.pos, m_pos);
    cv.force = detail::rebase_slice(host_cv.force, grid.force, m_force);
    cv.epsilon = detail::rebase_slice(host_cv.epsilon, grid.epsilon, m_epsilon);
    cv.orientation = detail::rebase_slice(
        host_cv.orientation, grid.orientation, m_orientation);
    return cv;
}

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<TurbineType, ActTrait>>>::
    copy_to_device()
{
    const auto& grid = m_data.grid();
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
        m_pos.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.force.begin(), grid.force.end(),
        m_force.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.epsilon.begin(), grid.epsilon.end(),
        m_epsilon.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.orientation.begin(),
        grid.orientation.end(), m_orientation.begin());

    // The component views stored in the metadata point to host memory, so
    // build equivalent views that point to the device copies above
    const auto& meta = m_data.meta();
    std::vector<ComponentView> blades_vec;
    blades_vec.reserve(meta.blades.size());
    for (const auto& cv : meta.blades) {
        blades_vec.push_back(to_device_view(cv));
    }
    std::vector<ComponentView> tower_vec(1, to_device_view(meta.tower));
    std::vector<ComponentView> hub_vec(1, to_device_view(meta.hub));

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, blades_vec.begin(), blades_vec.end(),
        m_blades.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tower_vec.begin(), tower_vec.end(),
        m_tower.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, hub_vec.begin(), hub_vec.end(),
        m_hub.begin());
}

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<TurbineType, ActTrait>>>::
operator()(const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    BL_PROFILE(
        "kynema-sgf::ActSrcOp<" + ActTrait::identifier() +
        ActSrcDisk::identifier() + ">");

    const auto& bx = mfi.tilebox();

    const auto bxa =
        kynema_sgf::utils::realbox_to_box(m_data.info().bound_box, geom);
    const auto& bxi = bx & bxa;
    if (bxi.isEmpty()) {
        return;
    }

    const auto& sarr = m_act_src(lev).array(mfi);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    auto& tdata = m_data.meta();
    const int nBlades = tdata.num_blades;
    const int nPB = tdata.num_pts_blade;
    const int nPT = tdata.num_pts_tower;
    const auto* blades = m_blades.data();
    const auto* tower = m_tower.data();
    const auto* hub = m_hub.data();
    // assume constant dr same for all blades
    const auto* host_pos = tdata.blades[0].pos.data();
    const amrex::Real dR = vs::mag_sqr(host_pos[0] - host_pos[1]);
    const amrex::Real dT =
        2.0_rt * std::numbers::pi_v<amrex::Real> / tdata.num_blades;
    // assume these will be copied to device by the lambda capture...
    const vs::Vector m_normal(tdata.rotor_frame.x());
    const vs::Vector m_origin(tdata.rot_center);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const vs::Vector cc{
            problo[0] + ((i + 0.5_rt) * dx[0]),
            problo[1] + ((j + 0.5_rt) * dx[1]),
            problo[2] + ((k + 0.5_rt) * dx[2]),
        };

        amrex::RealArray src_force = {0.0_rt};
        for (int ib = 0; ib < nBlades; ++ib) {
            // blade/disk contribution
            {
                const auto* pos = blades[ib].pos.data();
                const auto* force = blades[ib].force.data();
                const auto* eps = blades[ib].epsilon.data();

                for (int ip = 0; ip < nPB; ++ip) {
                    const auto R = utils::delta_pnts_cyl(
                                       m_origin, m_normal, m_origin, pos[ip])
                                       .x();
                    const auto dist_on_disk =
                        utils::delta_pnts_cyl(m_origin, m_normal, cc, pos[ip]);
                    const auto& pforce = force[ip];

                    const amrex::Real weight_R =
                        utils::linear_basis_1d(dist_on_disk.x(), dR);
                    const amrex::Real weight_T =
                        utils::linear_basis_1d(R * dist_on_disk.y(), dT * R);
                    const amrex::Real weight_N =
                        utils::gaussian1d(dist_on_disk.z(), eps[ip].x());
                    const auto projection_weight =
                        weight_R * weight_T * weight_N;

                    src_force[0] += projection_weight * pforce.x();
                    src_force[1] += projection_weight * pforce.y();
                    src_force[2] += projection_weight * pforce.z();
                }
            }
            // tower contribution
            {
                const auto* pos = tower[0].pos.data();
                const auto* force = tower[0].force.data();
                const auto* eps = tower[0].epsilon.data();
                const auto* tmat = tower[0].orientation.data();
                for (int ip = 0; ip < nPT; ++ip) {
                    const auto dist = cc - pos[ip];
                    const auto dist_local = tmat[ip] & dist;
                    const auto gauss_fac =
                        utils::gaussian3d(dist_local, eps[ip]);
                    const auto& pforce = force[ip];

                    src_force[0] += gauss_fac * pforce.x();
                    src_force[1] += gauss_fac * pforce.y();
                    src_force[2] += gauss_fac * pforce.z();
                }
            }
            // hub
            {
                const auto* pos = hub[0].pos.data();
                const auto* force = hub[0].force.data();
                const auto* eps = hub[0].epsilon.data();
                const auto* tmat = hub[0].orientation.data();
                const auto dist = cc - pos[0];
                const auto dist_local = tmat[0] & dist;
                const auto gauss_fac = utils::gaussian3d(dist_local, eps[0]);
                const auto& pforce = force[0];
                src_force[0] += gauss_fac * pforce.x();
                src_force[1] += gauss_fac * pforce.y();
                src_force[2] += gauss_fac * pforce.z();
            }
        }

        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}
} // namespace kynema_sgf::actuator::ops

#endif /* ACTSRCDISKOP_TURBINE_H */
