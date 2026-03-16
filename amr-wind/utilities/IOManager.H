#ifndef IOMANAGER_H
#define IOMANAGER_H

#include <string>
#include <unordered_map>
#include <set>

#include "AMReX_Vector.H"
#include "AMReX_BoxArray.H"
#include "AMReX_DistributionMapping.H"

namespace amr_wind {

class CFDSim;
class Field;
class IntField;
class DerivedQtyMgr;

/** Input/Output manager
 *  \ingroup utilities
 *
 *  This class contains infrastructure to output solution fields during the
 *  simulation as well as code that reads a checkpoint file during a restarted
 *  simulation. By default the \ref eqsys register the necessary fields during
 *  initialization and these fields are output into plot-file and checkpoint
 *  file formats based on output frequencies managed by SimTime. The user can
 *  request additional fields be output by setting appropriate parameters in the
 *  input file. The class also provides the ability to override output of the
 *  default fields and output a subset of those fields.
 */
class IOManager
{
public:
    explicit IOManager(CFDSim& /*sim*/);

    ~IOManager();

    IOManager(const IOManager&) = delete;
    IOManager& operator=(const IOManager&) = delete;

    //! Process user inputs and perform one-time initialization tasks
    void initialize_io();

    //! Write all user-requested fields to disk
    void write_plot_file();

    //! Write all necessary fields for restart
    void write_checkpoint_file(int start_level = 0, int end_level = -1);

    //! Read all necessary fields for a restart
    void read_checkpoint_fields(
        const std::string& restart_file,
        const amrex::Vector<amrex::BoxArray>& ba_chk,
        const amrex::Vector<amrex::DistributionMapping>& dm_chk,
        const amrex::IntVect& rep);

    //! Register a variable for output
    void register_output_var(const std::string& fname)
    {
        m_pltvars_default.insert(fname);
    }

    void register_output_int_var(const std::string& fname)
    {
        m_int_pltvars_default.insert(fname);
    }

    //! Register a variable for restart file
    void register_restart_var(const std::string& fname)
    {
        m_chkvars.insert(fname);
    }

    //! Convenience function to register a variable for both output and restart
    void register_io_var(const std::string& fname)
    {
        m_pltvars_default.insert(fname);
        m_chkvars.insert(fname);
    }

    //! Name of the checkpoint file to be read during a restart
    const std::string& restart_file() { return m_restart_file; }

    [[nodiscard]] bool is_restart() const { return !m_restart_file.empty(); }

    [[nodiscard]] const amrex::Vector<Field*>& plot_fields() const
    {
        return m_plt_fields;
    }

    [[nodiscard]] const amrex::Vector<Field*>& checkpoint_fields() const
    {
        return m_chk_fields;
    }

    //! Name of the post processing directory
    const std::string& post_processing_directory()
    {
        if (!amrex::UtilCreateDirectory(m_post_dir, 0755)) {
            amrex::CreateDirectoryFailed(m_post_dir);
        }
        return m_post_dir;
    }

private:
    void write_header(
        const std::string& /*chkname*/, int start_level, int end_level);

    void write_info_file(const std::string& /*path*/);

    CFDSim& m_sim;

    std::unique_ptr<DerivedQtyMgr> m_derived_mgr;

    //! Default output variables registered automatically in the code
    std::set<std::string> m_pltvars_default;

    //! Default output integer variables registered automatically in the code
    std::set<std::string> m_int_pltvars_default;

    //! Variables for output in checkpoint/restart files
    std::set<std::string> m_chkvars;

    //! Final list of fields to be output
    amrex::Vector<Field*> m_plt_fields;

    //! Final list of integer fields to be output
    amrex::Vector<IntField*> m_int_plt_fields;

    //! Final list of fields for restart
    amrex::Vector<Field*> m_chk_fields;

    //! Variable names (including components) for output
    amrex::Vector<std::string> m_plt_var_names;

    //! Prefix used for the plot file directories
    std::string m_plt_prefix{"plt"};

    //! Prefix used for the restart file directories
    std::string m_chk_prefix{"chk"};

    //! Name for the post_processing directory
    std::string m_post_dir{"post_processing"};

    //! Restart file name
    std::string m_restart_file;

    //! Total number of variables (including components) output to plot file
    int m_plt_num_comp{0};

    //! Flag indicating whether default fields should be output
    bool m_output_default_vars{true};

    //! Flag indicating whether we should allow missing restart fields
    bool m_allow_missing_restart_fields{true};

#ifdef AMR_WIND_USE_HDF5
    //! Flag indicating whether or not to output HDF5 plot files
    bool m_output_hdf5_plotfile{false};
#ifdef AMR_WIND_USE_HDF5_ZFP
    //! String for requesting a particular ZFP compression in HDF5 plot files
    std::string m_hdf5_compression{"None@0"};
#endif
#endif

    //! Number of plot and checkpoint data files per write
    int m_nfiles{256};
};

} // namespace amr_wind

#endif /* IOMANAGER_H */
