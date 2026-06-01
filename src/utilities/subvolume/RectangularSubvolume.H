#ifndef RECTANGULARSUBVOLUME_H
#define RECTANGULARSUBVOLUME_H

#include "src/utilities/subvolume/SubvolumeBase.H"

namespace kynema_sgf::subvolume {

/** Get subvolume defined by 3D rectangle
 *  \ingroup subvolume
 *
 *  Defines volume region with origin, number of points in each direction, and
 * target mesh resolution in each direction. Also has an optional parameter
 * regarding the data organization, which is the chunk size.
 */
class RectangularSubvolume
    : public SubvolumeBase::Register<RectangularSubvolume>
{
public:
    static std::string identifier() { return "Rectangular"; }

    explicit RectangularSubvolume(const CFDSim& /*sim*/);

    ~RectangularSubvolume() override;

    /** Read user inputs and initialize the subvolume object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Check if inputs have a valid correspondence in the mesh and save
    //! quantities for output steps
    void evaluate_inputs() override;

    //! Name of this subvolume object
    std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this subvolume object
    std::string subvolumetype() const override { return identifier(); }

    int lev() const override { return m_lev_for_sub; }

    amrex::BoxArray box_array() const override { return m_ba; }

protected:
    const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<int> m_npts_vec;
    amrex::Vector<amrex::Real> m_dx_vec;
    amrex::Vector<int> m_chunk_size_vec;

    int m_lev_for_sub{0};
    int m_verbose{0};

    bool m_chunk_size_present;

    amrex::BoxArray m_ba;

    std::string m_label;
};

} // namespace kynema_sgf::subvolume

#endif /* RECTANGULARSUBVOLUME_H */
