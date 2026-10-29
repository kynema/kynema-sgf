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

    // Device copies of the actuator grid data. The component views in the
    // turbine meta data hold host pointers and cannot be used in kernels.
    DeviceVecList m_pos;
    DeviceVecList m_force;
    DeviceVecList m_epsilon;
    DeviceTensorList m_orientation;

    //! Index of the first point of each component in the actuator grid
    int m_hub_start{0};
    int m_blade_start{0};
    int m_tower_start{0};

    void copy_to_device();

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
    const auto& meta = m_data.meta();
    m_pos.resize(grid.pos.size());
    m_force.resize(grid.force.size());
    m_epsilon.resize(grid.epsilon.size());
    m_orientation.resize(grid.orientation.size());

    // Locate the components within the actuator grid arrays
    const auto* gpos = grid.pos.data();
    m_hub_start = static_cast<int>(meta.hub.pos.data() - gpos);
    m_blade_start = static_cast<int>(meta.blades[0].pos.data() - gpos);
    for (int ib = 0; ib < meta.num_blades; ++ib) {
        AMREX_ALWAYS_ASSERT(
            meta.blades[ib].pos.data() ==
            gpos + m_blade_start + (ib * meta.num_pts_blade));
    }
    if (meta.num_pts_tower > 0) {
        m_tower_start = static_cast<int>(meta.tower.pos.data() - gpos);
    }
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
    const auto* pos = m_pos.data();
    const auto* force = m_force.data();
    const auto* eps = m_epsilon.data();
    const auto* tmat = m_orientation.data();
    const int hub_start = m_hub_start;
    const int blade_start = m_blade_start;
    const int tower_start = m_tower_start;
    // assume constant dr same for all blades
    const auto* host_pos = tdata.blades[0].pos.data();
    const amrex::Real dR = vs::mag(host_pos[0] - host_pos[1]);
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
        // blade/disk contribution
        for (int ib = 0; ib < nBlades; ++ib) {
            const int ioff = blade_start + (ib * nPB);
            for (int ip = ioff; ip < ioff + nPB; ++ip) {
                const auto R =
                    utils::delta_pnts_cyl(m_origin, m_normal, m_origin, pos[ip])
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
                const auto projection_weight = weight_R * weight_T * weight_N;

                src_force[0] += projection_weight * pforce.x();
                src_force[1] += projection_weight * pforce.y();
                src_force[2] += projection_weight * pforce.z();
            }
        }
        // tower contribution
        for (int ip = tower_start; ip < tower_start + nPT; ++ip) {
            const auto dist = cc - pos[ip];
            const auto dist_local = tmat[ip] & dist;
            const auto gauss_fac = utils::gaussian3d(dist_local, eps[ip]);
            const auto& pforce = force[ip];

            src_force[0] += gauss_fac * pforce.x();
            src_force[1] += gauss_fac * pforce.y();
            src_force[2] += gauss_fac * pforce.z();
        }
        // hub
        {
            const auto dist = cc - pos[hub_start];
            const auto dist_local = tmat[hub_start] & dist;
            const auto gauss_fac =
                utils::gaussian3d(dist_local, eps[hub_start]);
            const auto& pforce = force[hub_start];
            src_force[0] += gauss_fac * pforce.x();
            src_force[1] += gauss_fac * pforce.y();
            src_force[2] += gauss_fac * pforce.z();
        }

        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}
} // namespace kynema_sgf::actuator::ops

#endif /* ACTSRCDISKOP_TURBINE_H */
