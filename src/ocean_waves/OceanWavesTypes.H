#ifndef WAVEMODELTYPES_H
#define WAVEMODELTYPES_H

#include <utility>

#include "amr-wind/core/vs/vector_space.H"

#include "AMReX_Gpu.H"
#include "AMReX_RealBox.H"

namespace amr_wind {

class CFDSim;

namespace ocean_waves {

/** Abstract representation of an ocean wave type.
 *
 *  \ingroup immersed_boundary
 */
struct OceanWavesType
{};

using RealList = amrex::Vector<amrex::Real>;
using VecList = amrex::Vector<amr_wind::vs::Vector>;
using TensorList = amrex::Vector<amr_wind::vs::Tensor>;
using DeviceVecList = amrex::Gpu::DeviceVector<amr_wind::vs::Vector>;
using DeviceTensorList = amrex::Gpu::DeviceVector<amr_wind::vs::Tensor>;

/** Basic information that is common to all wave models
 *
 *  \ingroup ocean_waves
 */
struct OceanWavesInfo
{
    //! Unique string identifier for the ocean waves body (wing, turbine)
    std::string label;

    //! Bounding box used to determine processes where this ocean waves
    amrex::RealBox generation_region;

    amrex::RealBox absorption_refion;

    //! Unique integer identifier for the ocean waves
    int id{-1};

    OceanWavesInfo(std::string label_in, const int id_in)
        : label(std::move(label_in)), id(id_in)
    {}
};

/** Abstract representation of data holder for specific ocean waves types.
 *
 *  \ingroup ocean_waves
 *
 */
template <typename OceanWavesTrait>
class OceanWavesDataHolder
{
private:
    //! Instance of the CFD simulation environment
    CFDSim& m_sim;

    //! Basic information about this ocean waves component in relation to
    //! the simulation (usually IBInfo)
    typename OceanWavesTrait::InfoType m_info;

    //! Additional data necessary for a given wavem type
    typename OceanWavesTrait::MetaType m_meta;

public:
    /** Initialize the data structures
     *
     *  \param sim Instance of the CFD simulation environment
     *  \param label Unique string identifier for this immersed boundary
     * component \param id Unique global integer identifier for this immersed
     * boundary component
     */
    OceanWavesDataHolder(CFDSim& sim, const std::string& label, const int id)
        : m_sim(sim), m_info(label, id), m_meta()
    {}

    CFDSim& sim() { return m_sim; }
    [[nodiscard]] const CFDSim& sim() const { return m_sim; }

    typename OceanWavesTrait::InfoType& info() { return m_info; }
    [[nodiscard]] const typename OceanWavesTrait::InfoType& info() const
    {
        return m_info;
    }

    typename OceanWavesTrait::MetaType& meta() { return m_meta; }
    [[nodiscard]] const typename OceanWavesTrait::MetaType& meta() const
    {
        return m_meta;
    }
};

} // namespace ocean_waves
} // namespace amr_wind

#endif /* OCEANTYPES_H */
