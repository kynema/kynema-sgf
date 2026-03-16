#ifndef ACTUATOR_TYPES_H
#define ACTUATOR_TYPES_H

#include "amr-wind/core/Slice.H"
#include "amr-wind/core/vs/vector_space.H"

#include "AMReX_Gpu.H"
#include "AMReX_RealBox.H"

#include <set>
#include <utility>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class CFDSim;

namespace actuator {

/** Abstract representation of an actuator type.
 *
 *  \ingroup actuator
 */
struct ActuatorType
{};

/** Abstract representation of a source type.
 *
 *  \ingroup actuator
 */
struct ActSrcType
{};

/** Actuator line representation.
 *
 *  \ingroup actuator
 */
struct ActSrcLine : ActSrcType
{
    static std::string identifier() { return "Line"; }

    //! Flag indicating if this is a line type
    static constexpr bool is_line = true;

    //! Flag indicating if this is a disk type
    static constexpr bool is_disk = false;
};

/** Acutator disk representation.
 *
 *  \ingroup actuator
 */
struct ActSrcDisk : ActSrcType
{
    static std::string identifier() { return "Disk"; }

    static constexpr bool is_line = false;
    static constexpr bool is_disk = true;
};

using RealList = amrex::Vector<amrex::Real>;
using RealSlice = ::amr_wind::utils::Slice<amrex::Real>;
using VecList = amrex::Vector<amr_wind::vs::Vector>;
using VecSlice = ::amr_wind::utils::Slice<amr_wind::vs::Vector>;
using TensorList = amrex::Vector<amr_wind::vs::Tensor>;
using TensorSlice = ::amr_wind::utils::Slice<amr_wind::vs::Tensor>;
using DeviceVecList = amrex::Gpu::DeviceVector<amr_wind::vs::Vector>;
using DeviceTensorList = amrex::Gpu::DeviceVector<amr_wind::vs::Tensor>;

/** Actuator data at each node for an actuator component.
 *
 *  \ingroup actuator
 */
struct ActGrid
{
    //! Position vectors of the actuator forcing points on the grid
    VecList pos;

    //! Force vector at the forcing points
    VecList force;

    //! Gaussian smearing parameter
    VecList epsilon;

    //! Transformation matrix at the forcing points
    TensorList orientation;

    //! Position vectors for the points where velocity is sampled
    VecList vel_pos;

    //! Velocity vector at the sampled locations
    VecList vel;

    //! Density at the sampled locations
    RealList density;

    //! Switch for turning off distance components in Gaussian calculation
    vs::Vector dcoord_flags{1.0_rt, 1.0_rt, 1.0_rt};

    /** Helper method to resize the data arrays defined on the grid
     *
     *  \param num_force_pts Number of nodes that contain forcing data
     *  \param num_vel_pts Number of nodes where velocity field is sampled
     */
    void resize(int num_force_pts, int num_vel_pts)
    {
        pos.resize(num_force_pts);
        force.resize(num_force_pts);
        epsilon.resize(num_force_pts);
        orientation.resize(num_force_pts);
        vel_pos.resize(num_vel_pts);
        vel.resize(num_vel_pts);
        density.resize(num_vel_pts);
    }

    /** Convenience function to resize both force/velocity data to same
     * locations
     */
    void resize(int num_pts) { resize(num_pts, num_pts); }
};

/** Slice of actutator data for components (tower, blade, etc)
 *
 *  \ingroup actuator
 */
struct ComponentView
{
    VecSlice pos;
    VecSlice force;
    VecSlice epsilon;
    TensorSlice orientation;

    VecSlice vel_pos;
    VecSlice vel;
    VecSlice vel_rel;

    RealSlice chord;
};

/** Basic information that is common to all actuator types.
 *
 *  \ingroup actuator
 */
struct ActInfo
{
    //! Unique string identifier for the actuator body (wing, turbine)
    std::string label;

    //! Set of MPI ranks where this actuator source is active
    std::set<int> procs;

    //! Bounding box used to determine processes where this turbine has
    //! influence
    amrex::RealBox bound_box;

    //! Unique integer identifier for the turbine
    int id{-1};

    //! Root process where this turbine is active
    int root_proc{-1};

    //! Flag indicating whether this is root proc
    bool is_root_proc{false};

    //! Flag indicating whether this actuator component is active in the current
    //! MPI rank
    bool actuator_in_proc{false};

    //! Flag indicating whether this process requires velocities sampled at
    //! actuator point
    bool sample_vel_in_proc{false};

    ActInfo(std::string label_in, const int id_in)
        : label(std::move(label_in)), id(id_in)
    {}
};

/** Abstract representation of data holder for specific actuator types.
 *
 *  \ingroup actuator
 *
 */
template <typename ActTrait>
class ActDataHolder
{
private:
    //! Instance of the CFD simulation environment
    CFDSim& m_sim;

    //! Basic information about this actuator component in relation to the
    //! simulation (usually ActInfo)
    typename ActTrait::InfoType m_info;

    //! Nodal data on the actuator grid defined for this component (usually
    //! ActGrid)
    typename ActTrait::GridType m_grid;

    //! Additional data necessary for a given actuator type
    typename ActTrait::MetaType m_meta;

public:
    /** Initialize the data structures
     *
     *  \param sim Instance of the CFD simulation environment
     *  \param label Unique string identifier for this actuator component
     *  \param id Unique global integer identifier for this actuator component
     */
    ActDataHolder(CFDSim& sim, const std::string& label, const int id)
        : m_sim(sim), m_info(label, id), m_grid(), m_meta()
    {}

    CFDSim& sim() { return m_sim; }
    [[nodiscard]] const CFDSim& sim() const { return m_sim; }

    typename ActTrait::InfoType& info() { return m_info; }
    [[nodiscard]] const typename ActTrait::InfoType& info() const
    {
        return m_info;
    }

    typename ActTrait::GridType& grid() { return m_grid; }
    [[nodiscard]] const typename ActTrait::GridType& grid() const
    {
        return m_grid;
    }

    typename ActTrait::MetaType& meta() { return m_meta; }
    [[nodiscard]] const typename ActTrait::MetaType& meta() const
    {
        return m_meta;
    }
};

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATOR_TYPES_H */
