#ifndef ACTSRCDISKOP_H
#define ACTSRCDISKOP_H
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/wind_energy/actuator/disk/disk_spreading.H"

namespace amr_wind::actuator::ops {

template <typename ActTrait>
class ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<DiskType, ActTrait>>>
{
    using TraitType = ActTrait;

    template <typename T>
    using MyType = ActSrcOp<
        T,
        ActSrcDisk,
        std::enable_if_t<std::is_base_of_v<DiskType, T>>>;

    using SpreadType = SpreadingFunction<MyType<ActTrait>>;

    friend SpreadType;

protected:
    typename ActTrait::DataType& m_data;
    SpreadType m_spreading;
    Field& m_act_src;

    DeviceVecList m_pos;
    DeviceVecList m_force;

    void copy_to_device();

public:
    // cppcheck-suppress uninitMemberVar
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
    std::enable_if_t<std::is_base_of_v<DiskType, ActTrait>>>::
operator()(const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    BL_PROFILE(
        "amr-wind::ActSrcOp<" + ActTrait::identifier() +
        ActSrcDisk::identifier() + ">");
    m_spreading(*this, lev, mfi, geom);
}

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<DiskType, ActTrait>>>::initialize()
{
    const auto& grid = m_data.grid();
    m_pos.resize(grid.pos.size());
    m_force.resize(grid.force.size());
    m_spreading.initialize(m_data.meta().spreading_type);
}

template <typename ActTrait>
void ActSrcOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of_v<DiskType, ActTrait>>>::copy_to_device()
{
    const auto& grid = m_data.grid();

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.pos.begin(), grid.pos.end(),
        m_pos.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, grid.force.begin(), grid.force.end(),
        m_force.begin());
}

} // namespace amr_wind::actuator::ops

#endif /* ACTSRCDISKOP_H */
