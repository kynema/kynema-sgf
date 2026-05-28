#ifndef SUBVOLUMEBASE_H
#define SUBVOLUMEBASE_H

#include <AMReX_Geometry.H>
#include <AMReX_RealVect.H>
#include "src/core/Factory.H"

namespace kynema_sgf {

class CFDSim;

namespace subvolume {

/** Abstract representation of subvolume constraints to sample flow data
 *  \ingroup subvolume
 *
 *  This class defines the interface by which the positions of the subvolumes
 *  are defined and access from the Subvolume class. All subvolume definitions
 *  must be sub-classes of this base class.
 */
class SubvolumeBase : public Factory<SubvolumeBase, CFDSim&>
{
public:
    static std::string base_identifier() { return "SubvolumeBase"; }

    ~SubvolumeBase() override = default;

    //! Name used to refer to this sampler (e.g., myline1)
    [[nodiscard]] virtual std::string label() const = 0;
    virtual std::string& label() = 0;

    //! Class name of this sampler (e.g., Rectangular)
    [[nodiscard]] virtual std::string subvolumetype() const = 0;

    //! Read inputs
    virtual void initialize(const std::string& key) = 0;

    //! Check inputs and prepare for output
    virtual void evaluate_inputs() = 0;

    //! Run actions after initialize
    virtual void post_init_actions() { evaluate_inputs(); }

    //! Run actions after regrid
    virtual void post_regrid_actions() { evaluate_inputs(); }

    //! Level that the subvolume corresponds to
    [[nodiscard]] virtual int lev() const { return 0; }

    //! BoxArray that the subvolume corresponds to
    [[nodiscard]] virtual amrex::BoxArray box_array() const
    {
        return amrex::BoxArray(0);
    }

    //! Source BoxArray used to gather data before downsampling
    [[nodiscard]] virtual amrex::BoxArray source_box_array() const
    {
        return box_array();
    }

    //! Integer sampling stride in each direction from source cells
    [[nodiscard]] virtual amrex::IntVect sampling_stride() const
    {
        return amrex::IntVect{1, 1, 1};
    }

    //! Output geometry for the single-level plotfile
    [[nodiscard]] virtual amrex::Geometry output_geom() const
    {
        return amrex::Geometry{};
    }
};

} // namespace subvolume
} // namespace kynema_sgf

#endif /* SUBVOLUMEBASE_H */
