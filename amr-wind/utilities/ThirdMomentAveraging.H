#ifndef THIRDMOMENTAVERAGING_H
#define THIRDMOMENTAVERAGING_H

#include "AMReX_AmrCore.H"
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"

namespace amr_wind {

/** Compute the third moment with three fields
 *  \ingroup statistics
 *
 */
class ThirdMomentAveraging
{
public:
    ThirdMomentAveraging(
        FieldPlaneAveraging& pa1,
        FieldPlaneAveraging& pa2,
        FieldPlaneAveraging& pa3);

    ~ThirdMomentAveraging() = default;

    void operator()();

    /** evaluate third moment at specific location for both components */
    [[nodiscard]] amrex::Real line_average_interpolated(
        amrex::Real x, int comp1, int comp2, int comp3) const;
    /** evaluate third moment at specific location for each component */
    [[nodiscard]] amrex::Real
    line_average_interpolated(amrex::Real x, int comp) const;
    /** evaluate third moment at specific cell for both components */
    [[nodiscard]] amrex::Real
    line_average_cell(int ind, int comp1, int comp2, int comp3) const;
    /** evaluate third moment at specific cell for each component */
    [[nodiscard]] amrex::Real line_average_cell(int ind, int comp) const;

    const amrex::Vector<amrex::Real>& line_moment()
    {
        return m_third_moments_line;
    };
    void line_moment(int comp, amrex::Vector<amrex::Real>& l_vec);

    void output_line_average_ascii(
        const std::string& filename, int step, amrex::Real time);
    void output_line_average_ascii(int step, amrex::Real time);

    /** change precision of text file output */
    void set_precision(int p) { m_precision = p; };

private:
    int m_num_moments; /** outer product of components */
    amrex::Vector<amrex::Real>
        m_third_moments_line; /** line storage for fluctuations */

    int m_precision = 4;           /** precision for line plot text file */
    int m_last_updated_index = -1; /** keep track of the last time index that
                                      the operator was called */

    FieldPlaneAveraging& m_plane_average1;
    FieldPlaneAveraging& m_plane_average2;
    FieldPlaneAveraging& m_plane_average3;

public: // public for GPU
    /** fill line storage with averages */
    template <typename IndexSelector>
    void compute_average(
        const IndexSelector& idxOp,
        const amrex::MultiFab& mfab1,
        const amrex::MultiFab& mfab2,
        const amrex::MultiFab& mfab3);
};

} // namespace amr_wind

#endif /* ThirdMomentAveraging_H */
