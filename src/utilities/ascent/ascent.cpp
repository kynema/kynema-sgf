#include "ascent.H"

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/io_utils.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Conduit_Blueprint.H"

#include <ascent.hpp>

namespace amr_wind {
namespace ascent_int {

AscentPostProcess::AscentPostProcess(CFDSim& sim, const std::string& label)
    : m_sim(sim), m_label(label)
{}

AscentPostProcess::~AscentPostProcess() = default;

void AscentPostProcess::pre_init_actions() {}

void AscentPostProcess::initialize()
{
    BL_PROFILE("amr-wind::AscentPostProcess::initialize");

    amrex::Vector<std::string> field_names;

    {
        amrex::ParmParse pp("ascent");
        pp.getarr("fields", field_names);
        ioutils::assert_with_message(
            ioutils::all_distinct(field_names), "Duplicates in ascent.fields");
        populate_output_parameters(pp);
    }

    // Process field information
    auto& repo = m_sim.repo();

    for (const auto& fname : field_names) {
        if (!repo.field_exists(fname)) {
            amrex::Print() << "WARNING: Ascent: Non-existent field requested: "
                           << fname << std::endl;
            continue;
        }

        auto& fld = repo.get_field(fname);
        m_fields.emplace_back(&fld);
        ioutils::add_var_names(m_var_names, fld.name(), fld.num_comp());
    }
}

void AscentPostProcess::output_actions()
{
    BL_PROFILE("amr-wind::AscentPostProcess::output_actions");

    amrex::Vector<int> istep(
        m_sim.mesh().finestLevel() + 1, m_sim.time().time_index());

    int plt_num_comp = 0;
    for (auto* fld : m_fields) {
        plt_num_comp += fld->num_comp();
    }

    auto outfield = m_sim.repo().create_scratch_field(plt_num_comp);

    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        int icomp = 0;
        auto& mf = (*outfield)(lev);

        for (auto* fld : m_fields) {
            amrex::MultiFab::Copy(
                mf, (*fld)(lev), 0, icomp, fld->num_comp(), 0);
            icomp += fld->num_comp();
        }
    }

    const auto& mesh = m_sim.mesh();

    amrex::Print() << "Calling Ascent at time " << m_sim.time().new_time()
                   << std::endl;
    conduit::Node bp_mesh;
    amrex::MultiLevelToBlueprint(
        nlevels, outfield->vec_const_ptrs(), m_var_names, mesh.Geom(),
        m_sim.time().new_time(), istep, mesh.refRatio(), bp_mesh);

    ascent::Ascent ascent;
    conduit::Node open_opts;

#ifdef AMREX_USE_MPI
    open_opts["mpi_comm"] =
        MPI_Comm_c2f(amrex::ParallelDescriptor::Communicator());
#endif
    ascent.open(open_opts);
    conduit::Node verify_info;
    if (!conduit::blueprint::mesh::verify(bp_mesh, verify_info)) {
        ASCENT_INFO("Error: Mesh Blueprint Verify Failed!");
        verify_info.print();
    }

    conduit::Node actions;
    ascent.publish(bp_mesh);

    ascent.execute(actions);

    ascent.close();
}

void AscentPostProcess::post_regrid_actions()
{
    // nothing to do here
}

} // namespace ascent_int
} // namespace amr_wind
