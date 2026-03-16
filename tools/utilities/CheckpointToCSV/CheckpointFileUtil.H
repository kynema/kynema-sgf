#ifndef CHECKPOINTFILEUTIL_H_
#define CHECKPOINTFILEUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <string>
#include <memory>

// Adapted from amrex/Src/Base/AMReX_PlotFileDataImpl.H
// not inherited because constructor needed to be rewritten

class CheckpointFileDataImpl
{
public:
    explicit CheckpointFileDataImpl(
        std::string const& chkptfile_name,
        amrex::Vector<std::string> addl_field_names = {});

    [[nodiscard]] int spaceDim() const noexcept { return m_spacedim; }

    [[nodiscard]] amrex::Real time() const noexcept { return m_time; }

    [[nodiscard]] int finestLevel() const noexcept { return m_finest_level; }

    [[nodiscard]] const amrex::BoxArray& boxArray(int level) const noexcept
    {
        return m_ba[level];
    }

    [[nodiscard]] const amrex::DistributionMapping&
    DistributionMap(int level) const noexcept
    {
        return m_dmap[level];
    }

    void syncDistributionMap(CheckpointFileDataImpl const& src) noexcept;

    void
    syncDistributionMap(int level, CheckpointFileDataImpl const& src) noexcept;

    [[nodiscard]] amrex::Box probDomain(int level) const noexcept
    {
        return m_prob_domain[level];
    }
    [[nodiscard]] amrex::Array<amrex::Real, AMREX_SPACEDIM>
    probLo() const noexcept
    {
        return m_prob_lo;
    }
    [[nodiscard]] amrex::Array<amrex::Real, AMREX_SPACEDIM>
    probHi() const noexcept
    {
        return m_prob_hi;
    }
    [[nodiscard]] amrex::Array<amrex::Real, AMREX_SPACEDIM>
    cellSize(int level) const noexcept
    {
        return m_cell_size[level];
    }

    [[nodiscard]] const amrex::Vector<std::string>& varNames() const noexcept
    {
        return m_var_names;
    }

    [[nodiscard]] const amrex::Vector<std::string>&
    varNamesComponents() const noexcept
    {
        return m_var_names_full;
    }

    [[nodiscard]] int nComp() const noexcept { return m_ncomp; }
    [[nodiscard]] amrex::IntVect nGrowVect(int level) const noexcept
    {
        return m_ngrow[level];
    }

    amrex::MultiFab get(int level) noexcept;

private:
    std::string m_chkptfile_name;
    int m_ncomp{8};
    int m_nfields{4};
    amrex::Vector<std::string> m_var_names = {"density", "gp", "p", "velocity"};
    amrex::Vector<std::string> m_var_names_full = {
        "density", "gpx",       "gpy",       "gpz",
        "p",       "velocityx", "velocityy", "velocityz"};
    amrex::Vector<int> m_var_ncomp = {1, 3, 1, 3};
    const int m_spacedim{3};
    amrex::Real m_time;
    amrex::Real m_dt_restart;
    amrex::Real m_dt_nm1;
    amrex::Real m_dt_nm2;
    int m_finest_level, m_nlevels;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_prob_lo{
        {AMREX_D_DECL(0., 0., 0.)}};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_prob_hi{
        {AMREX_D_DECL(1., 1., 1.)}};
    amrex::Vector<int> m_ref_ratio;
    amrex::Vector<amrex::Box> m_prob_domain;
    int m_nstep;
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> m_cell_size;
    amrex::Vector<std::string> m_mf_name;
    amrex::Vector<std::unique_ptr<amrex::VisMF>> m_vismf;
    amrex::Vector<amrex::BoxArray> m_ba;
    amrex::Vector<amrex::DistributionMapping> m_dmap;
    amrex::Vector<amrex::IntVect> m_ngrow;
};

#endif
