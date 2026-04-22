#include "amr-wind/utilities/FieldPlaneAveraging.H"

#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

template <typename FType>
FPlaneAveraging<FType>::FPlaneAveraging(
    const FType& field_in,
    const amr_wind::SimTime& time,
    int axis_in,
    bool compute_deriv)
    : m_field(field_in)
    , m_time(time)
    , m_axis(axis_in)
    , m_comp_deriv(compute_deriv)
{
    AMREX_ALWAYS_ASSERT(m_axis >= 0 && m_axis < AMREX_SPACEDIM);
    auto geom = m_field.repo().mesh().Geom();
    // level=0 is default, could later make this an input.
    // Might only makes sense for fully covered levels
    m_xlo = geom[m_level].ProbLo(m_axis);
    m_dx = geom[m_level].CellSize(m_axis);

    m_ncomp = m_field.num_comp();

    const amrex::Box& domain = geom[m_level].Domain();
    const amrex::IntVect dom_lo(domain.loVect());
    const amrex::IntVect dom_hi(domain.hiVect());

    m_ncell_line = dom_hi[m_axis] - dom_lo[m_axis] + 1;

    // count number of cells in plane
    m_ncell_plane = 1;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (i != m_axis) {
            m_ncell_plane *= (dom_hi[i] - dom_lo[i] + 1);
        }
    }

    m_line_average.resize(static_cast<size_t>(m_ncell_line) * m_ncomp, 0.0_rt);
    if (m_comp_deriv) {
        m_line_deriv.resize(
            static_cast<size_t>(m_ncell_line) * m_ncomp, 0.0_rt);
    }
    m_line_xcentroid.resize(m_ncell_line);

    for (int i = 0; i < m_ncell_line; ++i) {
        m_line_xcentroid[i] = m_xlo + ((i + 0.5_rt) * m_dx);
    }
}

template <typename FType>
void FPlaneAveraging<FType>::output_line_average_ascii(
    const std::string& filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::FPlaneAveraging::output_line_average_ascii");

    if (step != m_last_updated_index) {
        operator()();
    }

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    std::ofstream outfile;
    outfile.precision(m_precision);

    if (step == 1) {
        // make new file
        outfile.open(filename.c_str(), std::ios_base::out);
        outfile << "#ncell,ncomp" << '\n';
        outfile << m_ncell_line << ", " << m_ncomp + 3 << '\n';
        outfile << "#step,time,z";
        for (int i = 0; i < m_ncomp; ++i) {
            outfile << ",<" + m_field.name() + std::to_string(i) + ">";
        }
        outfile << '\n';
    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
    }

    for (int i = 0; i < m_ncell_line; ++i) {
        outfile << step << ", " << std::scientific << time << ", "
                << m_line_xcentroid[i];
        for (int n = 0; n < m_ncomp; ++n) {
            outfile << ", " << std::scientific
                    << m_line_average[(m_ncomp * i) + n];
        }
        outfile << '\n';
    }
}

template <typename FType>
void FPlaneAveraging<FType>::output_line_average_ascii(
    int step, amrex::Real time)
{
    const std::string filename = "plane_average_" + m_field.name() + ".txt";
    output_line_average_ascii(filename, step, time);
}

template <typename FType>
amrex::Real
FPlaneAveraging<FType>::line_average_interpolated(amrex::Real x, int comp) const
{

    BL_PROFILE("amr-wind::PlaneAveraging::line_average_interpolated");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    amrex::Real c = 0.0_rt;
    int ind = 0;

    if (x > m_xlo + (0.5_rt * m_dx)) {
        ind = static_cast<int>(std::floor(((x - m_xlo) / m_dx) - 0.5_rt));
        const amrex::Real x1 = m_xlo + ((ind + 0.5_rt) * m_dx);
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0_rt;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);

    return (m_line_average[(m_ncomp * ind) + comp] * (1.0_rt - c)) +
           (m_line_average[(m_ncomp * (ind + 1)) + comp] * c);
}

template <typename FType>
void FPlaneAveraging<FType>::line_average(
    int comp, amrex::Vector<amrex::Real>& l_vec)
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_average");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    for (int i = 0; i < m_ncell_line; i++) {
        l_vec[i] = m_line_average[(m_ncomp * i) + comp];
    }
}

template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    return m_line_average[(m_ncomp * ind) + comp];
}

template <typename FType>
void FPlaneAveraging<FType>::compute_line_derivatives()
{
    BL_PROFILE("amr-wind::PlaneAveraging::compute_line_derivatives");
    for (int i = 0; i < m_ncell_line; ++i) {
        for (int n = 0; n < m_ncomp; ++n) {
            m_line_deriv[(m_ncomp * i) + n] =
                line_derivative_of_average_cell(i, n);
        }
    }
}

template <typename FType>
amrex::Real
FPlaneAveraging<FType>::line_derivative_of_average_cell(int ind, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_derivative_of_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    amrex::Real dudx;

    if (ind == 0) {
        dudx = (m_line_average[(m_ncomp * (ind + 1)) + comp] -
                m_line_average[(m_ncomp * ind) + comp]) /
               m_dx;
    } else if (ind == m_ncell_line - 1) {
        dudx = (m_line_average[(m_ncomp * (ind)) + comp] -
                m_line_average[(m_ncomp * (ind - 1)) + comp]) /
               m_dx;
    } else {
        dudx = 0.5_rt *
               (m_line_average[(m_ncomp * (ind + 1)) + comp] -
                m_line_average[(m_ncomp * (ind - 1)) + comp]) /
               m_dx;
    }

    return dudx;
}

template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_derivative_interpolated(
    amrex::Real x, int comp) const
{
    BL_PROFILE("amr-wind::PlaneAveraging::line_derivative_interpolated");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    amrex::Real c = 0.0_rt;
    int ind = 0;

    if (x > m_xlo + (0.5_rt * m_dx)) {
        ind = static_cast<int>(std::floor(((x - m_xlo) / m_dx) - 0.5_rt));
        const amrex::Real x1 = m_xlo + ((ind + 0.5_rt) * m_dx);
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0_rt;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);

    return (m_line_deriv[(m_ncomp * ind) + comp] * (1.0_rt - c)) +
           (m_line_deriv[(m_ncomp * (ind + 1)) + comp] * c);
}

template <typename FType>
void FPlaneAveraging<FType>::operator()()
{
    BL_PROFILE("amr-wind::FPlaneAveraging::operator");

    m_last_updated_index = m_time.time_index();

    std::ranges::fill(m_line_average, 0.0_rt);

    switch (m_axis) {
    case 0:
        compute_averages(XDir(), m_field(m_level));
        break;
    case 1:
        compute_averages(YDir(), m_field(m_level));
        break;
    case 2:
        compute_averages(ZDir(), m_field(m_level));
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }

    if (m_comp_deriv) {
        compute_line_derivatives();
    }
}

template <typename FType>
template <typename IndexSelector>
void FPlaneAveraging<FType>::compute_averages(
    const IndexSelector& idxOp, const amrex::MultiFab& mfab)
{
    BL_PROFILE("amr-wind::PlaneAveraging::compute_averages");

    const amrex::Real denom = 1.0_rt / (amrex::Real)m_ncell_plane;

    amrex::AsyncArray<amrex::Real> lavg(
        m_line_average.data(), m_line_average.size());

    amrex::Real* line_avg = lavg.data();
    const int captured_ncomp = m_ncomp;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        amrex::Box bx = mfi.tilebox();

        auto fab_arr = mfab.const_array(mfi);

        amrex::Box pbx =
            perpendicular_box<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

        amrex::ParallelFor(
            amrex::Gpu::KernelInfo().setReduction(true), pbx,
            [=] AMREX_GPU_DEVICE(
                int p_i, int p_j, int p_k, amrex::Gpu::Handler const& handler) {
                // Loop over the direction perpendicular to the plane.
                // This reduces the atomic pressure on the destination arrays.

                amrex::Box lbx = parallel_box<IndexSelector>(
                    bx, amrex::IntVect{p_i, p_j, p_k});

                for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                    for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                        for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0); ++i) {

                            const int ind = idxOp(i, j, k);

                            for (int n = 0; n < captured_ncomp; ++n) {
                                amrex::Gpu::deviceReduceSum(
                                    &line_avg[(captured_ncomp * ind) + n],
                                    fab_arr(i, j, k, n) * denom, handler);
                            }
                        }
                    }
                }
            });
    }

    lavg.copyToHost(m_line_average.data(), m_line_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_average.data(), m_line_average.size());

    // fixme later remove denom above and replace with this
    //    std::for_each(m_line_average.begin(), m_line_average.end(),
    //    [denom](amrex::Real &el){el *= denom; });
}

template class FPlaneAveraging<Field>;
template class FPlaneAveraging<ScratchField>;

// NOLINTBEGIN(clang-analyzer-security.ArrayBound)
VelPlaneAveraging::VelPlaneAveraging(CFDSim& sim, int axis_in)
    : FieldPlaneAveraging(
          sim.repo().get_field("velocity"), sim.time(), axis_in, true)
{
    m_line_hvelmag_average.resize(m_ncell_line, 0.0_rt);
    if (m_comp_deriv) {
        m_line_hvelmag_deriv.resize(m_ncell_line, 0.0_rt);
    }
}
// NOLINTEND(clang-analyzer-security.ArrayBound)

void VelPlaneAveraging::operator()()
{

    BL_PROFILE("amr-wind::VelPlaneAveraging::operator");

    FieldPlaneAveraging::operator()();

    std::ranges::fill(m_line_hvelmag_average, 0.0_rt);

    switch (m_axis) {
    case 0:
        compute_hvelmag_averages(XDir(), 1, 2, m_field(m_level));
        break;
    case 1:
        compute_hvelmag_averages(YDir(), 0, 2, m_field(m_level));
        break;
    case 2:
        compute_hvelmag_averages(ZDir(), 0, 1, m_field(m_level));
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }

    if (m_comp_deriv) {
        compute_line_hvelmag_derivatives();
    }
}

template <typename IndexSelector>
void VelPlaneAveraging::compute_hvelmag_averages(
    const IndexSelector& idx_op,
    const int h1_idx,
    const int h2_idx,
    const amrex::MultiFab& mfab)
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::compute_hvelmag_averages");

    const amrex::Real denom = 1.0_rt / (amrex::Real)m_ncell_plane;

    amrex::AsyncArray<amrex::Real> lavg(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::Real* line_avg = lavg.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        amrex::Box bx = mfi.tilebox();

        auto fab_arr = mfab.const_array(mfi);

        amrex::Box pbx =
            perpendicular_box<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

        amrex::ParallelFor(
            amrex::Gpu::KernelInfo().setReduction(true), pbx,
            [=] AMREX_GPU_DEVICE(
                int p_i, int p_j, int p_k, amrex::Gpu::Handler const& handler) {
                // Loop over the direction perpendicular to the plane.
                // This reduces the atomic pressure on the destination arrays.

                amrex::Box lbx = parallel_box<IndexSelector>(
                    bx, amrex::IntVect{p_i, p_j, p_k});

                for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                    for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                        for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0); ++i) {

                            const int ind = idx_op(i, j, k);
                            const amrex::Real hvelmag = std::sqrt(
                                (fab_arr(i, j, k, h1_idx) *
                                 fab_arr(i, j, k, h1_idx)) +
                                (fab_arr(i, j, k, h2_idx) *
                                 fab_arr(i, j, k, h2_idx)));
                            amrex::Gpu::deviceReduceSum(
                                &line_avg[ind], hvelmag * denom, handler);
                        }
                    }
                }
            });
    }

    lavg.copyToHost(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_hvelmag_average.data(),
        static_cast<int>(m_line_hvelmag_average.size()));
}

void VelPlaneAveraging::compute_line_hvelmag_derivatives()
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::compute_line_hvelmag_derivatives");
    for (int i = 0; i < m_ncell_line; ++i) {
        m_line_hvelmag_deriv[i] = line_hvelmag_derivative_of_average_cell(i);
    }
}

amrex::Real
VelPlaneAveraging::line_hvelmag_derivative_of_average_cell(int ind) const
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::line_derivative_of_average_cell");

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    amrex::Real dudx;

    if (ind == 0) {
        dudx =
            (m_line_hvelmag_average[(ind + 1)] - m_line_hvelmag_average[ind]) /
            m_dx;
    } else if (ind == m_ncell_line - 1) {
        dudx = (m_line_hvelmag_average[ind] - m_line_hvelmag_average[ind - 1]) /
               m_dx;
    } else {
        dudx = 0.5_rt *
               (m_line_hvelmag_average[ind + 1] -
                m_line_hvelmag_average[ind - 1]) /
               m_dx;
    }

    return dudx;
}

amrex::Real
VelPlaneAveraging::line_hvelmag_average_interpolated(amrex::Real x) const
{

    BL_PROFILE("amr-wind::PlaneAveraging::line_average_interpolated");

    amrex::Real c = 0.0_rt;
    int ind = 0;

    if (x > m_xlo + (0.5_rt * m_dx)) {
        ind = static_cast<int>(std::floor(((x - m_xlo) / m_dx) - 0.5_rt));
        const amrex::Real x1 = m_xlo + ((ind + 0.5_rt) * m_dx);
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0_rt;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);

    return (m_line_hvelmag_average[ind] * (1.0_rt - c)) +
           (m_line_hvelmag_average[ind + 1] * c);
}

amrex::Real VelPlaneAveraging::line_hvelmag_average_cell(int ind) const
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::line_hvelmag_average_cell");

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    return m_line_average[ind];
}

void VelPlaneAveraging::output_line_average_ascii(
    const std::string& filename, int step, amrex::Real time)
{
    BL_PROFILE("amr-wind::VelPlaneAveraging::output_line_average_ascii");

    if (step != m_last_updated_index) {
        operator()();
    }

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    std::ofstream outfile;
    outfile.precision(m_precision);

    if (step == 1) {
        // make new file
        outfile.open(filename.c_str(), std::ios_base::out);
        outfile << "#ncell,ncomp" << '\n';
        outfile << m_ncell_line << ", " << m_ncomp + 4 << '\n';
        outfile << "#step,time,z";
        for (int i = 0; i < m_ncomp; ++i) {
            outfile << ",<" + m_field.name() + std::to_string(i) + ">";
        }
        outfile << ", <hvelmag>";
        outfile << '\n';
    } else {
        // append file
        outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::app);
    }

    for (int i = 0; i < m_ncell_line; ++i) {
        outfile << step << ", " << std::scientific << time << ", "
                << m_line_xcentroid[i];
        for (int n = 0; n < m_ncomp; ++n) {
            outfile << ", " << std::scientific
                    << m_line_average[(m_ncomp * i) + n];
        }
        outfile << ", " << std::scientific << m_line_hvelmag_average[i];
        outfile << '\n';
    }
}

void VelPlaneAveraging::output_line_average_ascii(int step, amrex::Real time)
{
    const std::string filename = "plane_average_" + m_field.name() + ".txt";
    output_line_average_ascii(filename, step, time);
}

} // namespace amr_wind
