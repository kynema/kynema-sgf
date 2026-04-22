#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <chrono>
#include <ctime>
#include <fstream>

#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/DerivedQuantity.H"
#include "amr-wind/utilities/DerivedQtyDefs.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_MultiFabUtil.H"

#ifdef AMR_WIND_USE_HDF5
#include "AMReX_PlotFileUtilHDF5.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
#endif

namespace amr_wind {

IOManager::IOManager(CFDSim& sim)
    : m_sim(sim), m_derived_mgr(new DerivedQtyMgr(m_sim.repo()))
{}

IOManager::~IOManager() = default;

void IOManager::initialize_io()
{
    amrex::Vector<std::string> out_vars;
    amrex::Vector<std::string> out_int_vars;
    amrex::Vector<std::string> out_derived_vars;
    amrex::Vector<std::string> out_skip_vars;
    std::set<std::string> outputs;
    std::set<std::string> skip_outputs;
    std::set<std::string> int_outputs;

    amrex::ParmParse pp("io");
    pp.query("output_default_variables", m_output_default_vars);
    pp.query("plot_file", m_plt_prefix);
    pp.query("check_file", m_chk_prefix);
    pp.query("post_processing_directory", m_post_dir);
    pp.query("restart_file", m_restart_file);
    pp.query("allow_missing_restart_fields", m_allow_missing_restart_fields);
#ifdef AMR_WIND_USE_HDF5
    pp.query("output_hdf5_plotfile", m_output_hdf5_plotfile);
#ifdef AMR_WIND_USE_HDF5_ZFP
    pp.query("hdf5_compression", m_hdf5_compression);
#endif
#endif
    pp.query("nfiles", m_nfiles);

    // ParmParse requires us to read in a vector
    pp.queryarr("outputs", out_vars);
    pp.queryarr("int_outputs", out_int_vars);
    pp.queryarr("derived_outputs", out_derived_vars);
    pp.queryarr("skip_outputs", out_skip_vars);

    // We process the input vector to eliminate duplicates
    for (const auto& name : out_vars) {
        outputs.insert(name);
    }
    for (const auto& name : out_skip_vars) {
        skip_outputs.insert(name);
    }
    for (const auto& name : out_int_vars) {
        int_outputs.insert(name);
    }

    // If the user hasn't disabled default output variables, then we append them
    // to the list so that we don't end up with any duplicates (in case user
    // also added the variable explicitly in the input file)
    if (m_output_default_vars) {
        for (const auto& name : m_pltvars_default) {
            if (!skip_outputs.contains(name)) {
                outputs.insert(name);
            }
        }

        for (const auto& name : m_int_pltvars_default) {
            int_outputs.insert(name);
        }
    }

    amrex::Print() << "Initializing I/O manager" << '\n';

    // Process output variables information
    auto& repo = m_sim.repo();
    m_plt_num_comp = 0;
    for (const auto& fname : outputs) {
        if (repo.field_exists(fname)) {
            auto& fld = repo.get_field(fname);
            m_plt_num_comp += fld.num_comp();
            m_plt_fields.emplace_back(&fld);
            ioutils::add_var_names(m_plt_var_names, fld.name(), fld.num_comp());
        } else {
            amrex::Print() << "  Invalid output variable requested: " << fname
                           << '\n';
        }
    }

    for (const auto& fname : int_outputs) {
        if (repo.int_field_exists(fname)) {
            auto& fld = repo.get_int_field(fname);
            m_plt_num_comp += fld.num_comp();
            m_int_plt_fields.emplace_back(&fld);
            ioutils::add_var_names(m_plt_var_names, fld.name(), fld.num_comp());
        } else {
            amrex::Print() << "  Invalid output variable requested: " << fname
                           << '\n';
        }
    }

    if (!out_derived_vars.empty()) {
        m_derived_mgr->create(out_derived_vars);
        m_derived_mgr->filter(outputs);
        m_plt_num_comp += m_derived_mgr->num_comp();
        m_derived_mgr->var_names(m_plt_var_names);
    }

    for (const auto& fname : m_chkvars) {
        auto& fld = repo.get_field(fname);
        m_chk_fields.emplace_back(&fld);
    }

    amrex::VisMF::SetNOutFiles(m_nfiles);
}

void IOManager::write_plot_file()
{
    BL_PROFILE("amr-wind::IOManager::write_plot_file");

    amrex::Vector<int> istep(
        m_sim.mesh().finestLevel() + 1, m_sim.time().time_index());
    const int plt_comp = m_plt_num_comp;
    const int start_comp = m_plt_num_comp - m_derived_mgr->num_comp();
    auto outfield = m_sim.repo().create_scratch_field(plt_comp);
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        int icomp = 0;
        auto& mf = (*outfield)(lev);

        for (auto* fld : m_plt_fields) {
            amrex::MultiFab::Copy(
                mf, (*fld)(lev), 0, icomp, fld->num_comp(), 0);
            icomp += fld->num_comp();
        }

        for (auto* fld : m_int_plt_fields) {
            amrex::MultiFab::Copy(
                mf, amrex::ToMultiFab((*fld)(lev)), 0, icomp, fld->num_comp(),
                0);
            icomp += fld->num_comp();
        }
    }

    (*m_derived_mgr)(*outfield, start_comp);

    const std::string& plt_filename =
        amrex::Concatenate(m_plt_prefix, m_sim.time().time_index());
    const auto& mesh = m_sim.mesh();
    amrex::Print() << "Writing plot file       " << plt_filename << " at time "
                   << m_sim.time().new_time() << '\n';
#ifdef AMR_WIND_USE_HDF5
    if (m_output_hdf5_plotfile) {
        amrex::WriteMultiLevelPlotfileHDF5SingleDset(
            plt_filename, nlevels, outfield->vec_const_ptrs(), m_plt_var_names,
            mesh.Geom(), m_sim.time().new_time(), istep, mesh.refRatio()
#ifdef AMR_WIND_USE_HDF5_ZFP
                                                             ,
            m_hdf5_compression
#endif
        );
    } else {
#endif
        amrex::WriteMultiLevelPlotfile(
            plt_filename, nlevels, outfield->vec_const_ptrs(), m_plt_var_names,
            mesh.Geom(), m_sim.time().new_time(), istep, mesh.refRatio());
        write_info_file(plt_filename);
#ifdef AMR_WIND_USE_HDF5
    }
#endif
}

void IOManager::write_checkpoint_file(const int start_level, int end_level)
{
    BL_PROFILE("amr-wind::IOManager::write_checkpoint_file");
    const std::string level_prefix = "Level_";
    const std::string chkname =
        amrex::Concatenate(m_chk_prefix, m_sim.time().time_index());

    amrex::Print() << "Writing checkpoint file " << chkname << " at time "
                   << m_sim.time().new_time() << '\n';
    const auto& mesh = m_sim.mesh();
    // Modify end_level if need be
    end_level = (end_level == -1) ? mesh.finestLevel() : end_level;
    end_level = amrex::max(end_level, start_level);

    amrex::PreBuildDirectorHierarchy(
        chkname, level_prefix, end_level + 1 - start_level, true);
    write_header(chkname, start_level, end_level);
    write_info_file(chkname);

    for (int lev = start_level; lev < end_level + 1; ++lev) {
        for (auto* fld : m_chk_fields) {
            auto& field = *fld;
            amrex::VisMF::Write(
                field(lev),
                amrex::MultiFabFileFullPrefix(
                    lev - start_level, chkname, level_prefix, field.name()));
        }
    }
}

void IOManager::read_checkpoint_fields(
    const std::string& restart_file,
    const amrex::Vector<amrex::BoxArray>& ba_chk,
    const amrex::Vector<amrex::DistributionMapping>& dm_chk,
    const amrex::IntVect& rep)
{
    BL_PROFILE("amr-wind::IOManager::read_checkpoint_fields");

    // Track set of fields that might be missing at this level
    std::set<std::string> missing;
    const std::string level_prefix = "Level_";
    const int nlevels = m_sim.mesh().finestLevel() + 1;

    // always use the level 0 domain
    amrex::Box orig_domain(ba_chk[0].minimalBox());

    for (int lev = 0; lev < nlevels; ++lev) {
        for (auto* fld : m_chk_fields) {
            auto& field = *fld;
            const auto& fab_file = amrex::MultiFabFileFullPrefix(
                lev, restart_file, level_prefix, field.name());

            // Fields might be registered for checkpoint but might not be
            // necessary for actually performing the simulation. Check if the
            // field exists before attempting to read the restart field.
            if (!amrex::VisMF::Exist(fab_file)) {
                missing.insert(field.name());
                continue;
            }

            auto& mfab = field(lev);
            const auto& ba_fab = amrex::convert(ba_chk[lev], mfab.ixType());
            if (mfab.boxArray() == ba_fab &&
                mfab.DistributionMap() == dm_chk[lev]) {
                amrex::VisMF::Read(
                    field(lev),
                    amrex::MultiFabFileFullPrefix(
                        lev, restart_file, level_prefix, field.name()));
            } else {
                amrex::MultiFab tmp(
                    ba_fab, dm_chk[lev], mfab.nComp(), mfab.nGrowVect());
                amrex::VisMF::Read(
                    tmp, amrex::MultiFabFileFullPrefix(
                             lev, restart_file, level_prefix, field.name()));

                for (int k = 0; k < rep[2]; k++) {
                    for (int j = 0; j < rep[1]; j++) {
                        for (int i = 0; i < rep[0]; i++) {

                            amrex::IntVect shift_vec(
                                i * orig_domain.length(0),
                                j * orig_domain.length(1),
                                k * orig_domain.length(2));

                            // equivalent to 2^lev
                            shift_vec *= (1 << lev);

                            tmp.shift(shift_vec);
                            mfab.ParallelCopy(tmp);
                            tmp.shift(-shift_vec);
                        }
                    }
                }

                mfab.setBndry(0.0_rt);
            }
        }
    }

    // If fields were missing, print diagnostic message.
    if (!missing.empty()) {
        amrex::Print() << "\nWARNING: The following fields were missing in the "
                          "restart file for one or more levels. Please check "
                          "your restart file and inputs."
                       << '\n'
                       << "Missing checkpoint fields: " << '\n';
        for (const auto& ff : missing) {
            amrex::Print() << "  - " << ff << '\n';
        }
        amrex::Print() << '\n';
        if (!m_allow_missing_restart_fields) {
            amrex::Abort("Missing fields in restart file.");
        }
    }
}

void IOManager::write_header(
    const std::string& chkname, const int start_level, const int end_level)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string hdr_name(chkname + "/Header");
    amrex::VisMF::IO_Buffer io_buf(amrex::VisMF::IO_Buffer_Size);

    std::ofstream hdr;
    hdr.rdbuf()->pubsetbuf(io_buf.dataPtr(), io_buf.size());

    hdr.open(
        hdr_name.c_str(),
        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!hdr.good()) {
        amrex::FileOpenFailed(hdr_name);
    }

    hdr.precision(17);

    const auto& mesh = m_sim.mesh();
    const auto& time = m_sim.time();
    hdr << "Checkpoint version: 1\n"
        << end_level - start_level << "\n"
        << time.time_index() << "\n"
        << time.new_time() << "\n"
        << time.delta_t() << "\n"
        << time.delta_t_nm1() << "\n"
        << time.delta_t_nm2() << "\n";

    const auto geom = mesh.Geom(0);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        hdr << geom.ProbLo(i) << " ";
    }
    hdr << "\n";
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        hdr << geom.ProbHi(i) << " ";
    }
    hdr << "\n";

    for (int lev = start_level; lev < end_level + 1; ++lev) {
        mesh.boxArray(lev).writeOn(hdr);
        hdr << "\n";
    }

    hdr.close();
}

void IOManager::write_info_file(const std::string& path)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string dash_line = "\n" + std::string(78, '-') + "\n";
    const std::string fname(path + "/amr_wind_info");
    std::ofstream fh(fname.c_str(), std::ios::out);
    if (!fh.good()) {
        amrex::FileOpenFailed(fname);
    }

    amr_wind::io::print_banner(amrex::ParallelContext::CommunicatorSub(), fh);

    fh << dash_line << "Grid information: " << '\n';
    const auto& mesh = m_sim.mesh();
    for (int lev = 0; lev < m_sim.mesh().finestLevel() + 1; ++lev) {
        fh << "  Level: " << lev << "\n"
           << "    num. boxes = " << mesh.boxArray().size() << "\n"
           << "    maximum zones = ";

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            fh << mesh.Geom(lev).Domain().length(dir) << " ";
        }
        fh << "\n";
    }

    fh << dash_line << "Input file parameters: " << '\n';
    amrex::ParmParse::dumpTable(fh, true);
    fh.close();
}

} // namespace amr_wind
