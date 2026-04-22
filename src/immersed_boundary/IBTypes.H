#ifndef IBTYPES_H
#define IBTYPES_H

#include <utility>

#include "amr-wind/core/Slice.H"
#include "amr-wind/core/vs/vector_space.H"

#include "AMReX_Gpu.H"
#include "AMReX_RealBox.H"

namespace amr_wind {

class CFDSim;

namespace ib {

/** Abstract representation of an immersed boundary type.
 *
 *  \ingroup immersed_boundary
 */
struct IBType
{};

using RealList = amrex::Vector<amrex::Real>;
using RealSlice = ::amr_wind::utils::Slice<amrex::Real>;
using VecList = amrex::Vector<amr_wind::vs::Vector>;
using VecSlice = ::amr_wind::utils::Slice<amr_wind::vs::Vector>;
using TensorList = amrex::Vector<amr_wind::vs::Tensor>;
using TensorSlice = ::amr_wind::utils::Slice<amr_wind::vs::Tensor>;
using DeviceVecList = amrex::Gpu::DeviceVector<amr_wind::vs::Vector>;
using DeviceTensorList = amrex::Gpu::DeviceVector<amr_wind::vs::Tensor>;

/** Basic information that is common to all immersed boundary types.
 *
 *  \ingroup immersed_boundary
 */
struct IBInfo
{
    //! Unique string identifier for the immersed boundary body (wing, turbine)
    std::string label;

    //! Bounding box used to determine processes where this immersed boundary
    //! has influence
    amrex::RealBox bound_box;

    //! Unique integer identifier for the immersed boundary
    int id{-1};

    IBInfo(std::string label_in, const int id_in)
        : label(std::move(label_in)), id(id_in)
    {}
};

/** Abstract representation of data holder for specific immersed boundary types.
 *
 *  \ingroup immersed_boundary
 *
 */
template <typename IBTrait>
class IBDataHolder
{
private:
    //! Instance of the CFD simulation environment
    CFDSim& m_sim;

    //! Basic information about this immersed boundary component in relation to
    //! the simulation (usually IBInfo)
    typename IBTrait::InfoType m_info;

    //! Additional data necessary for a given actuator type
    typename IBTrait::MetaType m_meta;

public:
    /** Initialize the data structures
     *
     *  \param sim Instance of the CFD simulation environment
     *  \param label Unique string identifier for this immersed boundary
     * component \param id Unique global integer identifier for this immersed
     * boundary component
     */
    IBDataHolder(CFDSim& sim, const std::string& label, const int id)
        : m_sim(sim), m_info(label, id), m_meta()
    {}

    CFDSim& sim() { return m_sim; }
    [[nodiscard]] const CFDSim& sim() const { return m_sim; }

    typename IBTrait::InfoType& info() { return m_info; }
    [[nodiscard]] const typename IBTrait::InfoType& info() const
    {
        return m_info;
    }

    typename IBTrait::MetaType& meta() { return m_meta; }
    [[nodiscard]] const typename IBTrait::MetaType& meta() const
    {
        return m_meta;
    }
};

} // namespace ib
} // namespace amr_wind

#endif /* IBTYPES_H */
