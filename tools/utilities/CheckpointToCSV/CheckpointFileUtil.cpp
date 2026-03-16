#include "CheckpointFileUtil.H"
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_ParallelDescriptor.H>
#include <algorithm>
#include <fstream>
#include <iomanip>

namespace {
void GotoNextLine(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace

CheckpointFileDataImpl::CheckpointFileDataImpl(
    std::string const& chkptfile_name,
    amrex::Vector<std::string> addl_field_names)
    : m_chkptfile_name(chkptfile_name)
{
    // Add field names from input
    const int naf = addl_field_names.size();
    for (int n = 0; n < naf; ++n) {
        m_var_names.push_back(addl_field_names[n]);
        m_var_names_full.push_back(addl_field_names[n]);
        m_var_ncomp.push_back(1);
    }
    // Total number of components, assume addl fields are one component
    m_ncomp += naf;
    // Total number of fields
    m_nfields += naf;

    // Header
    std::string File(chkptfile_name + "/Header");
    amrex::Vector<char> fileCharPtr;
    amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::istringstream is(
        std::string(fileCharPtr.dataPtr()), std::istringstream::in);

    std::string line, word;

    // !! -- Start reading from checkpoint file -- !!
    // Title line
    std::getline(is, line);

    // Finest level
    is >> m_finest_level;
    GotoNextLine(is);

    // Step count
    is >> m_nstep;
    GotoNextLine(is);

    // Current time
    is >> m_time;
    GotoNextLine(is);

    // Time step size
    is >> m_dt_restart;
    GotoNextLine(is);

    is >> m_dt_nm1;
    GotoNextLine(is);

    is >> m_dt_nm2;
    GotoNextLine(is);

    // Low coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_prob_lo[i++] = static_cast<amrex::Real>(std::stod(word));
        }
    }

    // High coordinates of domain bounding box
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_prob_hi[i++] = static_cast<amrex::Real>(std::stod(word));
        }
    }

    // Read in box array to get domain box
    m_nlevels = m_finest_level + 1;
    m_ba.resize(m_nlevels);
    m_dmap.resize(m_nlevels);
    m_ngrow.resize(m_nlevels);
    m_prob_domain.resize(m_nlevels);
    m_cell_size.resize(
        m_nlevels,
        amrex::Array<amrex::Real, AMREX_SPACEDIM>{{AMREX_D_DECL(1., 1., 1.)}});
    for (int lev = 0; lev <= m_finest_level; ++lev) {
        // read in level 'lev' BoxArray from Header
        m_ba[lev].readFrom(is);
        GotoNextLine(is);
        // get minimal box from box array for prob domain
        m_prob_domain[lev] = m_ba[lev].minimalBox();
        // make distribution map
        m_dmap[lev].define(m_ba[lev], amrex::ParallelDescriptor::NProcs());
        // ngrow is set to 0, don't know how to find it properly
        m_ngrow[lev] = {0, 0, 0};
    }

    // Calculate cell size for lowest level
    for (int idim = 0; idim < m_spacedim; ++idim) {
        m_cell_size[0][idim] = (m_prob_hi[idim] - m_prob_lo[idim]) /
                               (m_prob_domain[0].bigEnd(idim) -
                                m_prob_domain[0].smallEnd(idim) + 1);
    }
    // Do the higher levels using refratio, assume 2
    for (int lev = 1; lev <= m_finest_level; ++lev) {
        for (int idim = 0; idim < m_spacedim; ++idim) {
            m_cell_size[lev][idim] = m_cell_size[lev - 1][idim];
        }
    }

    AMREX_ASSERT(m_nlevels > 0 && m_nlevels <= 1000);
}

void CheckpointFileDataImpl::syncDistributionMap(
    CheckpointFileDataImpl const& src) noexcept
{
    int nlevs_min = std::min(m_nlevels, src.m_nlevels);
    for (int ilev = 0; ilev < nlevs_min; ++ilev) {
        syncDistributionMap(ilev, src);
    }
}

void CheckpointFileDataImpl::syncDistributionMap(
    int level, CheckpointFileDataImpl const& src) noexcept
{
    if (level <= src.finestLevel() &&
        m_dmap[level].size() == src.DistributionMap(level).size()) {
        m_dmap[level] = src.DistributionMap(level);
    }
}

amrex::MultiFab CheckpointFileDataImpl::get(int level) noexcept
{
    const std::string level_prefix{"Level_"};
    amrex::MultiFab mf(m_ba[level], m_dmap[level], m_ncomp, m_ngrow[level]);

    // Do checkpoint reading, which is a field at a time
    for (int lev = 0; lev < m_nlevels; ++lev) {
        int icomp = 0;
        for (int nf = 0; nf < m_nfields; ++nf) {
            amrex::MultiFab tmp_mfab(
                m_ba[level], m_dmap[level], m_var_ncomp[nf], m_ngrow[level]);

            // Read current field into temporary fab
            amrex::VisMF::Read(
                tmp_mfab,
                amrex::MultiFabFileFullPrefix(
                    lev, m_chkptfile_name, level_prefix, m_var_names[nf]));

            // Copy from tmp fab to bigger fab
            amrex::MultiFab::Copy(
                mf, tmp_mfab, 0, icomp, m_var_ncomp[nf], m_ngrow[lev]);
            // Advance index among total components
            icomp += m_var_ncomp[nf];
        }
    }
    return mf;
}
