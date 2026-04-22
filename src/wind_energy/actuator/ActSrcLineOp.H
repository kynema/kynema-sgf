#ifndef ACTSRCLINEOP_H
#define ACTSRCLINEOP_H

#include "AMReX_REAL.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/core/FieldRepo.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops {

template <typename ActTrait>
class ActSrcOp<ActTrait, ActSrcLine>
{
private:
    typename ActTrait::DataType& m_data;
    Field& m_act_src;

    DeviceVecList m_pos;
    DeviceVecList m_pos_old;
    DeviceVecList m_force;
    DeviceVecList m_epsilon;
    DeviceTensorList m_orientation;

    bool m_init_old{false};

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
void ActSrcOp<ActTrait, ActSrcLine>::initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_pos_old.resize(grid.pos.size());
    m_force.resize(grid.force.size());
    m_epsilon.resize(grid.epsilon.size());
    m_orientation.resize(grid.orientation.size());
}

template <typename ActTrait>
void ActSrcOp<ActTrait, ActSrcLine>::copy_to_device()
{
    const auto& grid = m_data.grid();

    // Populate old positions before new positions are updated
    if (m_init_old) {
        amrex::Gpu::copy(
            amrex::Gpu::deviceToDevice, m_pos.begin(), m_pos.end(),
            m_pos_old.begin());
    }

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

    if (!m_init_old) {
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
            m_pos_old.begin());
        m_init_old = true;
    }
}

template <typename ActTrait>
void ActSrcOp<ActTrait, ActSrcLine>::operator()(
    const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    BL_PROFILE(
        "amr-wind::ActSrcOp<" + ActTrait::identifier() +
        ActSrcLine::identifier() + ">");

    const auto& bx = mfi.tilebox();

    const auto bxa =
        amr_wind::utils::realbox_to_box(m_data.info().bound_box, geom);
    const auto& bxi = bx & bxa;
    if (bxi.isEmpty()) {
        return;
    }

    const auto& sarr = m_act_src(lev).array(mfi);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    const int npts = m_pos.size();
    const auto* pos = m_pos.data();
    const auto* opos = m_pos_old.data();
    const auto* force = m_force.data();
    const auto* eps = m_epsilon.data();
    const auto* tmat = m_orientation.data();

    const auto dcoord_flags = m_data.grid().dcoord_flags;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const vs::Vector cc{
            problo[0] + ((i + 0.5_rt) * dx[0]),
            problo[1] + ((j + 0.5_rt) * dx[1]),
            problo[2] + ((k + 0.5_rt) * dx[2]),
        };

        amrex::RealArray src_force = {0.0_rt};
        for (int ip = 0; ip < npts; ++ip) {
            // Put force at n+1/2 location for Godunov
            constexpr amrex::Real wt = 0.5_rt;
            const auto pos_ip = wt * pos[ip] + (1.0_rt - wt) * opos[ip];
            const auto dist = cc - pos_ip;
            // Convert to local (chord, span, thickness) coords
            const auto dist_local_3D = tmat[ip] & dist;
            // In local coords, zero disabled distances (e.g., for 2D)
            const auto dist_local = dist_local_3D * dcoord_flags;
            const auto gauss_fac = utils::gaussian3d(dist_local, eps[ip]);
            const auto& pforce = force[ip];

            src_force[0] += gauss_fac * pforce.x();
            src_force[1] += gauss_fac * pforce.y();
            src_force[2] += gauss_fac * pforce.z();
        }

        sarr(i, j, k, 0) += src_force[0];
        sarr(i, j, k, 1) += src_force[1];
        sarr(i, j, k, 2) += src_force[2];
    });
}
} // namespace amr_wind::actuator::ops

#endif /* ACTSRCLINEOP_H */
