#include "src/utilities/FieldPlaneAveraging.H"
#include "src/utilities/constants.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

template <typename FType>
FPlaneAveraging<FType>::FPlaneAveraging(
    const FType& field_in,
    const kynema_sgf::SimTime& time,
    int axis_in,
    int max_level,
    bool compute_deriv)
    : m_field(field_in)
    , m_time(time)
    , m_axis(axis_in)
    , m_max_level(max_level)
    , m_comp_deriv(compute_deriv)
{
    AMREX_ALWAYS_ASSERT(m_axis >= 0 && m_axis < AMREX_SPACEDIM);
    const auto& mesh = m_field.repo().mesh();
    auto geom = mesh.Geom();

    // beginning and end of line, for now assuming line is the length of the
    // entire domain
    m_xlo = geom[0].ProbLo(m_axis);
    m_xhi = geom[0].ProbHi(m_axis);

    int finestLevel = mesh.maxLevel();
    if (m_max_level >= 0) {
        finestLevel = std::min(m_max_level, finestLevel);
    }
    const amrex::IntVect dom_lo_vec(geom[finestLevel].Domain().smallEnd());
    const amrex::IntVect dom_hi_vec(geom[finestLevel].Domain().bigEnd());
    const int dom_lo = dom_lo_vec[m_axis];
    const int dom_hi = dom_hi_vec[m_axis];

    AMREX_ALWAYS_ASSERT(dom_lo == 0);
    int dom_hi2 = geom[0].Domain().bigEnd()[m_axis] + 1;
    for (int i = 0; i < finestLevel; ++i) {
        dom_hi2 *= mesh.refRatio(i)[m_axis];
    }
    AMREX_ALWAYS_ASSERT(dom_hi + 1 == dom_hi2);

    // TODO: make an input maybe?
    m_ncell_line = dom_hi - dom_lo + 1;

    m_dx = (m_xhi - m_xlo) / static_cast<amrex::Real>(m_ncell_line);

    m_ncomp = m_field.num_comp();

    // count number of cells in plane
    m_ncell_plane = 1;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (i != m_axis) {
            m_ncell_plane *= (dom_hi_vec[i] - dom_lo_vec[i] + 1);
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
void FPlaneAveraging<FType>::convert_x_to_ind(
    amrex::Real x, int& ind, amrex::Real& c) const
{
    c = 0.0_rt;
    ind = 0;

    if (x > m_xlo + 0.5_rt * m_dx) {
        ind = static_cast<int>(std::floor((x - m_xlo) / m_dx - 0.5_rt));
        const amrex::Real x1 = m_xlo + (ind + 0.5_rt) * m_dx;
        c = (x - x1) / m_dx;
    }

    if (ind + 1 >= m_ncell_line) {
        ind = m_ncell_line - 2;
        c = 1.0_rt;
    }

    AMREX_ALWAYS_ASSERT(ind >= 0 && ind + 1 < m_ncell_line);
}

template <typename FType>
void FPlaneAveraging<FType>::output_line_average_ascii(
    const std::string& filename, int step, amrex::Real time)
{
    BL_PROFILE("kynema-sgf::FPlaneAveraging::output_line_average_ascii");

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

    BL_PROFILE("kynema-sgf::PlaneAveraging::line_average_interpolated");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return (m_line_average[(m_ncomp * ind) + comp] * (1.0_rt - c)) +
           (m_line_average[(m_ncomp * (ind + 1)) + comp] * c);
}

template <typename FType>
void FPlaneAveraging<FType>::line_average(
    int comp, amrex::Vector<amrex::Real>& l_vec)
{
    BL_PROFILE("kynema-sgf::PlaneAveraging::line_average");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    for (int i = 0; i < m_ncell_line; i++) {
        l_vec[i] = m_line_average[(m_ncomp * i) + comp];
    }
}

template <typename FType>
amrex::Real FPlaneAveraging<FType>::line_average_cell(int ind, int comp) const
{
    BL_PROFILE("kynema-sgf::PlaneAveraging::line_average_cell");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);
    AMREX_ALWAYS_ASSERT(ind >= 0 && ind < m_ncell_line);

    return m_line_average[(m_ncomp * ind) + comp];
}

template <typename FType>
void FPlaneAveraging<FType>::operator()()
{
    BL_PROFILE("kynema-sgf::FPlaneAveraging::operator");

    m_last_updated_index = m_time.time_index();

    std::fill(m_line_average.begin(), m_line_average.end(), 0.0_rt);

    switch (m_axis) {
    case 0:
        compute_averages(XDir());
        break;
    case 1:
        compute_averages(YDir());
        break;
    case 2:
        compute_averages(ZDir());
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
void FPlaneAveraging<FType>::compute_averages(const IndexSelector& idxOp)
{
    BL_PROFILE("kynema-sgf::PlaneAveraging::compute_averages");

    amrex::AsyncArray<amrex::Real> lavg(
        m_line_average.data(), m_line_average.size());

    amrex::Real* line_avg = lavg.data();
    const int num_comps = m_ncomp;
    const int num_cells = m_ncell_line;
    const amrex::Real line_dx = m_dx;
    const amrex::Real xlo = m_xlo;

    auto g0 = m_field.repo().mesh().Geom(0);
    amrex::Real lateral_area = 1.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != m_axis) {
            lateral_area *= (g0.ProbHi(d) - g0.ProbLo(d));
        }
    }
    const amrex::Real denom = 1.0_rt / (lateral_area * m_dx);

    const bool periodic_dir = g0.periodicity().isPeriodic(m_axis);

    const auto& mesh = m_field.repo().mesh();
    int finestLevel = mesh.finestLevel();
    if (m_max_level >= 0) {
        finestLevel = std::min(m_max_level, finestLevel);
    }
    const auto dir = m_axis;
    const bool no_ghost = (m_field.num_grow()[dir] == 0);

    for (int lev = 0; lev <= finestLevel; ++lev) {

        const auto& geom = mesh.Geom(lev);
        const amrex::Real dx = geom.CellSize()[dir];
        const amrex::Real dy = geom.CellSize()[idxOp.odir1];
        const amrex::Real dz = geom.CellSize()[idxOp.odir2];

        const amrex::Real problo_x = geom.ProbLo(dir);
        const amrex::Real probhi_x = geom.ProbHi(dir);

        amrex::iMultiFab level_mask;
        if (lev < finestLevel) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& mfab = m_field(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            amrex::Box bx = mfi.tilebox();
            amrex::Box vbx = mfi.validbox();

            auto fab_arr = mfab.const_array(mfi);
            auto mask_arr = level_mask.const_array(mfi);

            amrex::Box pbx =
                perpendicular_box<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

            amrex::ParallelFor(
                amrex::Gpu::KernelInfo().setReduction(true), pbx,
                [=] AMREX_GPU_DEVICE(
                    int p_i, int p_j, int p_k,
                    amrex::Gpu::Handler const& handler) noexcept {
                    // Loop over the direction perpendicular to the plane.
                    // This reduces the atomic pressure on the destination
                    // arrays.

                    amrex::Box lbx = parallel_box<IndexSelector>(
                        bx, amrex::IntVect{p_i, p_j, p_k});

                    for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                        for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                            for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0);
                                 ++i) {

                                // cell coordinates
                                const amrex::Real cell_xlo =
                                    xlo + idxOp(i, j, k) * dx;
                                const amrex::Real cell_xhi = cell_xlo + dx;

                                // line indices
                                const int line_ind_lo = amrex::min(
                                    amrex::max(
                                        static_cast<int>(
                                            (cell_xlo - xlo) / line_dx),
                                        0),
                                    num_cells - 1);
                                const int line_ind_hi = amrex::min(
                                    amrex::max(
                                        static_cast<int>(
                                            (cell_xhi -
                                             kynema_sgf::constants::TIGHT_TOL -
                                             xlo) /
                                            line_dx),
                                        0),
                                    num_cells - 1);

                                AMREX_ASSERT(line_ind_lo >= 0);
                                AMREX_ASSERT(line_ind_hi >= 0);
                                AMREX_ASSERT(line_ind_lo < num_cells);
                                AMREX_ASSERT(line_ind_hi < num_cells);

                                for (int ind = line_ind_lo; ind <= line_ind_hi;
                                     ++ind) {

                                    // line coordinates
                                    const amrex::Real line_xlo =
                                        xlo + ind * line_dx;
                                    const amrex::Real line_xhi =
                                        line_xlo + line_dx;

                                    amrex::Real deltax;

                                    if (line_xlo <= cell_xlo) {
                                        deltax = line_xhi - cell_xlo;
                                    } else if (line_xhi >= cell_xhi) {
                                        deltax = cell_xhi - line_xlo;
                                    } else {
                                        deltax = line_dx;
                                    }
                                    deltax = amrex::min(deltax, dx);
                                    const amrex::Real vol = deltax * dy * dz;

                                    // Calculate location of target
                                    const auto x_targ =
                                        0.5_rt * (line_xlo + line_xhi);
                                    // Calculate location of cell center
                                    const amrex::IntVect iv{i, j, k};
                                    const auto idx = iv[dir];
                                    const auto x_cell =
                                        problo_x +
                                        (static_cast<amrex::Real>(idx) +
                                         0.5_rt) *
                                            dx;
                                    // Get location of neighboring cell centers
                                    auto x_up = x_cell + dx;
                                    auto x_down = x_cell - dx;
                                    // Bound locations by domain limits
                                    if (!periodic_dir) {
                                        x_up = amrex::min(probhi_x, x_up);
                                        x_down = amrex::max(problo_x, x_down);
                                    }
                                    // Pick indices of closest neighbor
                                    // Bound indices in case of no ghost cells
                                    auto iv_nb = iv;
                                    auto x_nb = x_cell;
                                    if (std::abs(x_up - x_targ) <
                                        std::abs(x_down - x_targ)) {
                                        x_nb = x_up;
                                        iv_nb[dir] += 1;
                                        if (no_ghost) {
                                            iv_nb[dir] = amrex::min<int>(
                                                iv_nb[dir], vbx.bigEnd(dir));
                                        }
                                    } else {
                                        x_nb = x_down;
                                        iv_nb[dir] -= 1;
                                        if (no_ghost) {
                                            iv_nb[dir] = amrex::max<int>(
                                                iv_nb[dir], vbx.smallEnd(dir));
                                        }
                                    }
                                    // Interpolate to target location using
                                    // closest neighbor
                                    // (will do nothing if already at cell
                                    // center)
                                    for (int n = 0; n < num_comps; ++n) {
                                        const auto f_interp =
                                            fab_arr(iv, n) +
                                            (fab_arr(iv_nb, n) -
                                             fab_arr(iv, n)) *
                                                ((x_targ - x_cell) /
                                                 (x_nb - x_cell));
                                        amrex::Gpu::deviceReduceSum(
                                            &line_avg[num_comps * ind + n],
                                            mask_arr(i, j, k) * f_interp * vol *
                                                denom,
                                            handler);
                                    }
                                }
                            }
                        }
                    }
                });
        }
    }

    lavg.copyToHost(m_line_average.data(), m_line_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_average.data(), m_line_average.size());
}

template <typename FType>
void FPlaneAveraging<FType>::compute_line_derivatives()
{
    BL_PROFILE("kynema-sgf::PlaneAveraging::compute_line_derivatives");
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
    BL_PROFILE("kynema-sgf::PlaneAveraging::line_derivative_of_average_cell");

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
    BL_PROFILE("kynema-sgf::PlaneAveraging::line_derivative_interpolated");

    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < m_ncomp);

    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return (m_line_deriv[(m_ncomp * ind) + comp] * (1.0_rt - c)) +
           (m_line_deriv[(m_ncomp * (ind + 1)) + comp] * c);
}

template class FPlaneAveraging<Field>;
template class FPlaneAveraging<ScratchField>;

VelPlaneAveraging::VelPlaneAveraging(CFDSim& sim, int axis_in, int max_level)
    : FieldPlaneAveraging(
          sim.repo().get_field("velocity"),
          sim.time(),
          axis_in,
          max_level,
          true)
{
    m_line_hvelmag_average.resize(m_ncell_line, 0.0_rt);
    m_line_Su_average.resize(m_ncell_line, 0.0_rt);
    m_line_Sv_average.resize(m_ncell_line, 0.0_rt);
}
// NOLINTEND(clang-analyzer-security.ArrayBound)

void VelPlaneAveraging::operator()()
{

    BL_PROFILE("kynema-sgf::VelPlaneAveraging::operator");

    // velocity averages
    FieldPlaneAveraging::operator()();

    std::fill(
        m_line_hvelmag_average.begin(), m_line_hvelmag_average.end(), 0.0_rt);
    std::fill(m_line_Su_average.begin(), m_line_Su_average.end(), 0.0_rt);
    std::fill(m_line_Sv_average.begin(), m_line_Sv_average.end(), 0.0_rt);

    switch (m_axis) {
    case 0:
        compute_hvelmag_averages(XDir());
        break;
    case 1:
        compute_hvelmag_averages(YDir());
        break;
    case 2:
        compute_hvelmag_averages(ZDir());
        break;
    default:
        amrex::Abort("axis must be equal to 0, 1, or 2");
        break;
    }
}

template <typename IndexSelector>
void VelPlaneAveraging::compute_hvelmag_averages(const IndexSelector& idxOp)
{
    BL_PROFILE("kynema-sgf::VelPlaneAveraging::compute_hvelmag_averages");

    amrex::AsyncArray<amrex::Real> lavg_vm(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    amrex::AsyncArray<amrex::Real> lavg_Su(
        m_line_Su_average.data(), m_line_Su_average.size());
    amrex::AsyncArray<amrex::Real> lavg_Sv(
        m_line_Sv_average.data(), m_line_Sv_average.size());

    amrex::Real* line_avg_vm = lavg_vm.data();
    amrex::Real* line_avg_Su = lavg_Su.data();
    amrex::Real* line_avg_Sv = lavg_Sv.data();

    const int num_cells = m_ncell_line;
    const amrex::Real line_dx = m_dx;
    const amrex::Real xlo = m_xlo;

    const auto& mesh = m_field.repo().mesh();
    auto g0 = mesh.Geom(0);
    amrex::Real lateral_area = 1.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != m_axis) {
            lateral_area *= (g0.ProbHi(d) - g0.ProbLo(d));
        }
    }
    const amrex::Real denom = 1.0_rt / (lateral_area * m_dx);

    int finestLevel = mesh.finestLevel();
    if (m_max_level >= 0) {
        finestLevel = std::min(m_max_level, finestLevel);
    }

    for (int lev = 0; lev <= finestLevel; ++lev) {

        const auto& geom = mesh.Geom(lev);
        const amrex::Real dx = geom.CellSize()[m_axis];
        const amrex::Real dy = geom.CellSize()[idxOp.odir1];
        const amrex::Real dz = geom.CellSize()[idxOp.odir2];

        amrex::iMultiFab level_mask;
        if (lev < finestLevel) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        const auto& mfab = m_field(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            amrex::Box bx = mfi.tilebox();

            auto fab_arr = mfab.const_array(mfi);
            auto mask_arr = level_mask.const_array(mfi);

            amrex::Box pbx =
                perpendicular_box<IndexSelector>(bx, amrex::IntVect{0, 0, 0});

            amrex::ParallelFor(
                amrex::Gpu::KernelInfo().setReduction(true), pbx,
                [=] AMREX_GPU_DEVICE(
                    int p_i, int p_j, int p_k,
                    amrex::Gpu::Handler const& handler) noexcept {
                    // Loop over the direction perpendicular to the plane.
                    // This reduces the atomic pressure on the destination
                    // arrays.

                    amrex::Box lbx = parallel_box<IndexSelector>(
                        bx, amrex::IntVect{p_i, p_j, p_k});

                    for (int k = lbx.smallEnd(2); k <= lbx.bigEnd(2); ++k) {
                        for (int j = lbx.smallEnd(1); j <= lbx.bigEnd(1); ++j) {
                            for (int i = lbx.smallEnd(0); i <= lbx.bigEnd(0);
                                 ++i) {

                                // cell coordinates
                                const amrex::Real cell_xlo =
                                    xlo + idxOp(i, j, k) * dx;
                                const amrex::Real cell_xhi = cell_xlo + dx;

                                // line indices
                                const int line_ind_lo = amrex::min(
                                    amrex::max(
                                        static_cast<int>(
                                            (cell_xlo - xlo) / line_dx),
                                        0),
                                    num_cells - 1);
                                const int line_ind_hi = amrex::min(
                                    amrex::max(
                                        static_cast<int>(
                                            (cell_xhi - xlo) / line_dx),
                                        0),
                                    num_cells - 1);

                                AMREX_ASSERT(line_ind_lo >= 0);
                                AMREX_ASSERT(line_ind_hi >= 0);
                                AMREX_ASSERT(line_ind_lo < num_cells);
                                AMREX_ASSERT(line_ind_hi < num_cells);

                                for (int ind = line_ind_lo; ind <= line_ind_hi;
                                     ++ind) {

                                    // line coordinates
                                    const amrex::Real line_xlo =
                                        xlo + ind * line_dx;
                                    const amrex::Real line_xhi =
                                        line_xlo + line_dx;

                                    amrex::Real deltax;

                                    if (line_xlo <= cell_xlo) {
                                        deltax = line_xhi - cell_xlo;
                                    } else if (line_xhi >= cell_xhi) {
                                        deltax = cell_xhi - line_xlo;
                                    } else {
                                        deltax = line_dx;
                                    }

                                    deltax = amrex::min(deltax, dx);
                                    const amrex::Real vol =
                                        mask_arr(i, j, k) * deltax * dy * dz;

                                    const amrex::Real hvelmag = std::sqrt(
                                        fab_arr(i, j, k, idxOp.odir1) *
                                            fab_arr(i, j, k, idxOp.odir1) +
                                        fab_arr(i, j, k, idxOp.odir2) *
                                            fab_arr(i, j, k, idxOp.odir2));
                                    const amrex::Real Su =
                                        hvelmag * fab_arr(i, j, k, idxOp.odir1);
                                    const amrex::Real Sv =
                                        hvelmag * fab_arr(i, j, k, idxOp.odir2);

                                    amrex::Gpu::deviceReduceSum(
                                        &line_avg_vm[ind],
                                        hvelmag * vol * denom, handler);
                                    amrex::Gpu::deviceReduceSum(
                                        &line_avg_Su[ind], Su * vol * denom,
                                        handler);
                                    amrex::Gpu::deviceReduceSum(
                                        &line_avg_Sv[ind], Sv * vol * denom,
                                        handler);
                                }
                            }
                        }
                    }
                });
        }
    }

    lavg_vm.copyToHost(
        m_line_hvelmag_average.data(), m_line_hvelmag_average.size());
    lavg_Su.copyToHost(m_line_Su_average.data(), m_line_Su_average.size());
    lavg_Sv.copyToHost(m_line_Sv_average.data(), m_line_Sv_average.size());
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_hvelmag_average.data(),
        static_cast<int>(m_line_hvelmag_average.size()));
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_Su_average.data(), static_cast<int>(m_line_Su_average.size()));
    amrex::ParallelDescriptor::ReduceRealSum(
        m_line_Sv_average.data(), static_cast<int>(m_line_Sv_average.size()));
}

amrex::Real
VelPlaneAveraging::line_hvelmag_average_interpolated(amrex::Real x) const
{
    BL_PROFILE(
        "kynema-sgf::VelPlaneAveraging::line_hvelmag_average_interpolated");
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return (m_line_hvelmag_average[ind] * (1.0_rt - c)) +
           (m_line_hvelmag_average[ind + 1] * c);
}

amrex::Real VelPlaneAveraging::line_su_average_interpolated(amrex::Real x) const
{
    BL_PROFILE("kynema-sgf::VelPlaneAveraging::line_su_average_interpolated");
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_Su_average[ind] * (1.0_rt - c) +
           m_line_Su_average[ind + 1] * c;
}

amrex::Real VelPlaneAveraging::line_sv_average_interpolated(amrex::Real x) const
{
    BL_PROFILE("kynema-sgf::VelPlaneAveraging::line_sv_average_interpolated");
    int ind;
    amrex::Real c;
    convert_x_to_ind(x, ind, c);

    return m_line_Sv_average[ind] * (1.0_rt - c) +
           m_line_Sv_average[ind + 1] * c;
}

} // namespace kynema_sgf
