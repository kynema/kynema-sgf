#ifndef SECONDMOMENTAVERAGING_H
#define SECONDMOMENTAVERAGING_H

#include <AMReX_AmrCore.H>
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"

namespace amr_wind {

/** Compute second moments for two variables
 *  \ingroup statistics
 *
 *  Computes second-moment statistics for any two variables. For example, if the
 *  velocity field is provided, then it computes the following quantities:
 *
 *  \f$\langle u'u'\rangle, \langle u'v'\rangle, \langle u'w'\rangle, \ldots\f$
 *
 *  Alternately, using velocity (vector) and a scalar (e.g., temperature
 *  \f$\theta\f$), would compute:
 *
 *  \f$\langle u'\theta ' \rangle, \langle v' \theta ' \rangle, \langle w'
 * \theta ' \rangle\f$
 */
class SecondMomentAveraging
{
public:
    SecondMomentAveraging(FieldPlaneAveraging& pa1, FieldPlaneAveraging& pa2);

    ~SecondMomentAveraging() = default;

    void operator()();

    /** evaluate second moment at specific location for both components */
    [[nodiscard]] amrex::Real
    line_average_interpolated(amrex::Real x, int comp1, int comp2) const;
    /** evaluate second moment at specific location for each component */
    [[nodiscard]] amrex::Real
    line_average_interpolated(amrex::Real x, int comp) const;
    /** evaluate second moment at specific cell for both components */
    [[nodiscard]] amrex::Real
    line_average_cell(int ind, int comp1, int comp2) const;
    /** evaluate second moment at specific cell for each component */
    [[nodiscard]] amrex::Real line_average_cell(int ind, int comp) const;

    const amrex::Vector<amrex::Real>& line_moment()
    {
        return m_second_moments_line;
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
        m_second_moments_line; /** line storage for fluctuations */

    int m_precision = 4;           /** precision for line plot text file */
    int m_last_updated_index = -1; /** keep track of the last time index that
                                      the operator was called */

    FieldPlaneAveraging& m_plane_average1;
    FieldPlaneAveraging& m_plane_average2;

public: // public for GPU
    /** fill line storage with averages */
    template <typename IndexSelector>
    void compute_average(
        const IndexSelector& idxOp,
        const amrex::MultiFab& mfab1,
        const amrex::MultiFab& mfab2);
};

} // namespace amr_wind

#endif /* SecondMomentAveraging_H */
