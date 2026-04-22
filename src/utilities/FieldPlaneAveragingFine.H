#ifndef FIELDPLANEAVERAGINGFINE_H
#define FIELDPLANEAVERAGINGFINE_H

#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/SimTime.H"

/**
 *  \defgroup statistics Field statistics
 *  Field statistics
 *
 *  This group contains utilities for performing turbulence averaging and
 *  outputting statistics during wind simulations.
 *
 *  \ingroup utilities
 */

namespace amr_wind {

/** Output average of a field on planes normal to a given direction
 *  \ingroup statistics we_abl
 *
 *  The user can choose a direction (x, y, or z), the default value is the
 *  z-direction. The field is then averaged on planes at the cell-centers at
 *  level 0 (coarsest level) along the specified direction.
 */

template <typename FType>
class FPlaneAveragingFine
{
public:
    /**
     *  \param field_in [in] Field to be averaged
     *  \param time [in] Time instance to determine output frequencies
     *  \param axis_in [in] Direction along which planes are computed
     *  \param compute_deriv [in] Should the derivative of the averages be
     * computed
     */
    FPlaneAveragingFine(
        const FType& field_in,
        const amr_wind::SimTime& time,
        int axis_in,
        bool compute_deriv = false);

    virtual ~FPlaneAveragingFine() = default;

    virtual void operator()();

    void convert_x_to_ind(amrex::Real x, int& ind, amrex::Real& c) const;

    /** evaluate line average at specific location for any average component */
    [[nodiscard]] amrex::Real
    line_average_interpolated(amrex::Real x, int comp) const;
    /** evaluate line average at specific cell for any average component */
    [[nodiscard]] amrex::Real line_average_cell(int ind, int comp) const;

    /** evaluate line average derivative at specific location for any average
        component */
    [[nodiscard]] amrex::Real
    line_derivative_interpolated(amrex::Real x, int comp) const;
    /** evaluate derivative of a line average at specific cell for any component
     */
    [[nodiscard]] amrex::Real
    line_derivative_of_average_cell(int ind, int comp) const;

    void output_line_average_ascii(
        const std::string& filename, int step, amrex::Real time);
    void output_line_average_ascii(int step, amrex::Real time);

    /** change precision of text file output */
    void set_precision(int p) { m_precision = p; };

    [[nodiscard]] amrex::Real dx() const { return m_dx; };
    [[nodiscard]] amrex::Real xlo() const { return m_xlo; };

    [[nodiscard]] int axis() const { return m_axis; };
    [[nodiscard]] int ncomp() const { return m_ncomp; };
    [[nodiscard]] int ncell_line() const { return m_ncell_line; };
    [[nodiscard]] int last_updated_index() const
    {
        return m_last_updated_index;
    };

    [[nodiscard]] const amrex::Vector<amrex::Real>& line_average() const
    {
        return m_line_average;
    };
    [[nodiscard]] const amrex::Vector<amrex::Real>& line_deriv() const
    {
        return m_line_deriv;
    };
    void line_average(int comp, amrex::Vector<amrex::Real>& l_vec);
    [[nodiscard]] const amrex::Vector<amrex::Real>& line_centroids() const
    {
        return m_line_xcentroid;
    };

    [[nodiscard]] const FType& field() const { return m_field; };

protected:
    int m_ncomp; /** number of average components */

    amrex::Vector<amrex::Real>
        m_line_average; /** line storage for the average velocity and tracer
                           variables */
    //! line storage for the derivative of average field
    amrex::Vector<amrex::Real> m_line_deriv;

    amrex::Vector<amrex::Real> m_line_xcentroid; /** line storage for centroids
                                                    of each cell along a line*/

    amrex::Real m_dx;  /** mesh spacing in axis direction*/
    amrex::Real m_xlo; /** bottom of line */
    amrex::Real m_xhi; /** top of line */

    int m_ncell_line; /** number of cells along line */

    int m_precision = 6;           /** precision for line plot text file */
    int m_last_updated_index = -1; /** keep track of the last time index that
                                      the operator was called */

    const FType& m_field;
    const SimTime& m_time;
    const int m_axis;
    const bool m_comp_deriv;

public: // public for GPU
    /** fill line storage with averages */
    template <typename IndexSelector>
    void compute_averages(const IndexSelector& idxOp);

    /** fill derivatives of line storage */
    void compute_line_derivatives();
};

using FieldPlaneAveragingFine = FPlaneAveragingFine<Field>;
using ScratchFieldPlaneAveragingFine = FPlaneAveragingFine<ScratchField>;

/** Specific application of FieldPlaneAveraging to the velocity field
 *  \ingroup statistics we_abl
 */
class VelPlaneAveragingFine : public FieldPlaneAveragingFine
{
public:
    VelPlaneAveragingFine(CFDSim& sim, int axis_in);

    ~VelPlaneAveragingFine() override = default;

    void operator()() override;

private:
    amrex::Vector<amrex::Real>
        m_line_hvelmag_average; /** line storage for the average horizontal
                                   velocity magnitude */
    amrex::Vector<amrex::Real>
        m_line_Su_average; /** line storage for the average horizontal
                                   velocity magnitude time x-velocity */

    amrex::Vector<amrex::Real>
        m_line_Sv_average; /** line storage for the average horizontal
                                   velocity magnitude time y-velocity */

public: // public for GPU
    /** fill line storage with horizontal velocity magnitude averages */
    template <typename IndexSelector>
    void compute_hvelmag_averages(const IndexSelector& idxOp);

    /** evaluate line haverage at specific location for horizontal velocity
     * magnitude */
    [[nodiscard]] amrex::Real
    line_hvelmag_average_interpolated(amrex::Real x) const;
    /** evaluate line haverage at specific location for horizontal velocity
     * magnitude times x-velocity */
    [[nodiscard]] amrex::Real line_su_average_interpolated(amrex::Real x) const;
    /** evaluate line haverage at specific location for horizontal velocity
     * magnitude times y-velocity */
    [[nodiscard]] amrex::Real line_sv_average_interpolated(amrex::Real x) const;
};

} // namespace amr_wind

#endif /* FieldPlaneAveragingFine_H */
