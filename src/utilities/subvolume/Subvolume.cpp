#include <memory>
#include <utility>

#include "src/utilities/subvolume/Subvolume.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/IOManager.H"
#include "AMReX_MultiFabUtil.H"

#include "AMReX_ParmParse.H"

namespace kynema_sgf::subvolume {

Subvolume::Subvolume(CFDSim& sim, std::string label)
    : m_sim(sim)
    , m_derived_mgr(new DerivedQtyMgr(m_sim.repo()))
    , m_label(std::move(label))
{}

Subvolume::~Subvolume() = default;

void Subvolume::initialize()
{
    BL_PROFILE("amr-wind::Subvolume::initialize");

    // Labels for the different sampler types
    amrex::Vector<std::string> labels;
    // Fields to be sampled - requested by user
    amrex::Vector<std::string> field_names;
    // Int Fields to be sampled - requested by user
    amrex::Vector<std::string> int_field_names;
    // Derived fields to be sampled - requested by user
    amrex::Vector<std::string> derived_field_names;

    {
        amrex::ParmParse pp(m_label);
        populate_output_parameters(pp);
        pp.getarr("labels", labels);
        ioutils::assert_with_message(
            ioutils::all_distinct(labels),
            "Duplicates in " + m_label + ".labels");
        pp.getarr("fields", field_names);
        ioutils::assert_with_message(
            ioutils::all_distinct(field_names),
            "Duplicates in " + m_label + ".fields");
        pp.queryarr("int_fields", int_field_names);
        ioutils::assert_with_message(
            ioutils::all_distinct(int_field_names),
            "Duplicates in " + m_label + ".int_fields");
        pp.queryarr("derived_fields", derived_field_names);
        // Label can be renamed for the sake of output
        m_label_out = m_label;
        pp.query("output_rename", m_label_out);
        if (m_label_out != m_label) {
            amrex::Print()
                << "WARNING: Subvolume top-level label " << m_label
                << " renamed to " << m_label_out
                << " for output. User must ensure that bottom-level labels are "
                   "not redundant across different subvolume instances.\n";
        }
    }

    // Process field information
    m_ncomp = 0;
    auto& repo = m_sim.repo();
    for (const auto& fname : field_names) {
        if (!repo.field_exists(fname)) {
            amrex::Print()
                << "WARNING: Subvolume: Non-existent field requested: " << fname
                << ". This is a mistake or the requested field is a int field "
                   "or a derived field and should be added to the "
                   "int_fields/derived_fields parameter.\n";
            continue;
        }

        auto& fld = repo.get_field(fname);
        m_ncomp += fld.num_comp();
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    // Process field information
    m_nicomp = 0;
    for (const auto& fname : int_field_names) {
        if (!repo.int_field_exists(fname)) {
            amrex::Print()
                << "WARNING: Subvolume: Non-existent int_field requested: "
                << fname
                << ". This is a mistake or the requested int_field is a "
                   "derived field and should be added to the derived_fields "
                   "parameter\n";
            continue;
        }

        auto& fld = repo.get_int_field(fname);
        m_nicomp += fld.num_comp();
        m_int_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }

    // Process derived field information
    if (!derived_field_names.empty()) {
        m_derived_mgr->create(derived_field_names);
        m_derived_mgr->filter(field_names);
        m_ndcomp = m_derived_mgr->num_comp();
        m_derived_mgr->var_names(m_var_names);
    }

    // Load Subvolumes, only type is Rectangular right now
    for (const auto& lbl : labels) {
        const std::string key = m_label + "." + lbl;
        amrex::ParmParse pp1(key);
        std::string stype = "Rectangular";

        pp1.query("type", stype);
        auto obj = SubvolumeBase::create(stype, m_sim);
        obj->label() = lbl;
        obj->initialize(key);
        obj->post_init_actions();

        m_subvolumes.emplace_back(std::move(obj));
    }
}

void Subvolume::post_regrid_actions()
{
    BL_PROFILE("amr-wind::Subvolume::post_regrid_actions");

    for (const auto& obj : m_subvolumes) {
        obj->post_regrid_actions();
    }
}

void Subvolume::output_actions()
{
    BL_PROFILE("amr-wind::Subvolume::output_actions");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string name(post_dir + "/" + m_label_out);

    const auto time = m_sim.time().new_time();
    const auto itime = m_sim.time().time_index();
    const auto n_out = m_ncomp + m_nicomp + m_ndcomp;

    // Create scratch field for derived quantities
    auto scr_ptr = m_sim.repo().create_scratch_field(amrex::max(1, m_ndcomp));
    // Populate it
    if (m_ndcomp > 0) {
        (*m_derived_mgr)(*scr_ptr, 0);
    }

    // Current approach
    // - Populate scratch field only for sake of derived quantities
    // - Loop through subvolumes
    // -- Copy fields of interest to mf of entire level
    // -- Then copy from this mf to output mf_sv

    // Alternative way of doing this:
    // - Loop through subvolumes to get levels of all subvolumes
    // - Loop through these valid levels
    // -- Copy fields of interest to scratch field, now used for all outputs
    // - Loop through subvolumes again
    // -- Copy from scratch field to output mf_sv
    // ** would prevent a lot of redundant copies if many subvolumes are on the
    //    same level

    for (const auto& sv : m_subvolumes) {

        const auto lev = sv->lev();
        const auto ba = sv->box_array();

        amrex::DistributionMapping dm(ba);
        amrex::MultiFab mf_sv(ba, dm, n_out, 0);

        amrex::MultiFab mf_all_samelev(
            m_sim.mesh().boxArray(lev), m_sim.mesh().DistributionMap(lev),
            n_out, 0);

        int icomp = 0;
        for (auto* fld : m_fields) {
            amrex::MultiFab::Copy(
                mf_all_samelev, (*fld)(lev), 0, icomp, fld->num_comp(), 0);
            icomp += fld->num_comp();
        }

        for (auto* fld : m_int_fields) {
            amrex::MultiFab::Copy(
                mf_all_samelev, amrex::ToMultiFab((*fld)(lev)), 0, icomp,
                fld->num_comp(), 0);
            icomp += fld->num_comp();
        }

        // Contribute derived quantities to mf_all_samelev
        if (m_ndcomp > 0) {
            amrex::MultiFab::Copy(
                mf_all_samelev, (*scr_ptr)(lev), 0, icomp, m_ndcomp, 0);
        }

        mf_sv.ParallelCopy(mf_all_samelev, 0, 0, n_out, 0, 0);

        std::string sv_label = name + "_" + sv->label();
        std::string subvol_filename =
            amrex::Concatenate(sv_label, m_sim.time().time_index(), 5);

        amrex::Print() << "Writing subvolume into " << subvol_filename << "\n";
        WriteSingleLevelPlotfile(
            subvol_filename, mf_sv, m_var_names, m_sim.mesh().Geom(lev), time,
            itime);
    }
}

} // namespace kynema_sgf::subvolume
