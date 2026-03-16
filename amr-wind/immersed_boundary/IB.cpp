#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/IBModel.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/io_utils.H"
#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ib {

IB::IB(CFDSim& sim)
    : m_sim(sim)
    , m_ib_levelset(sim.repo().declare_field("ib_levelset", 1, 1, 1))
    , m_ib_normal(sim.repo().declare_field("ib_normal", AMREX_SPACEDIM, 1, 1))
{
    m_ib_levelset.set_default_fillpatch_bc(sim.time());
    m_ib_normal.set_default_fillpatch_bc(sim.time());
}

IB::~IB() = default;

void IB::pre_init_actions()
{
    BL_PROFILE("amr-wind::ib::IB::pre_init_actions");
    amrex::ParmParse pp(identifier());

    amrex::Vector<std::string> labels;
    pp.getarr("labels", labels);
    ioutils::assert_with_message(
        ioutils::all_distinct(labels),
        "Duplicates in " + identifier() + ".labels");

    const int n_ibs = static_cast<int>(labels.size());

    for (int i = 0; i < n_ibs; ++i) {
        const std::string& tname = labels[i];
        const std::string& prefix = identifier() + "." + tname;
        amrex::ParmParse pp1(prefix);

        std::string type;
        pp.query("type", type);
        pp1.query("type", type);
        AMREX_ALWAYS_ASSERT(!type.empty());

        auto obj = ImmersedBoundaryModel::create(type, m_sim, tname, i);

        const std::string default_prefix = identifier() + "." + type;
        ::amr_wind::utils::MultiParser inp(default_prefix, prefix);

        obj->read_inputs(inp);
        m_ibs.emplace_back(std::move(obj));
    }
}

void IB::post_init_actions()
{
    BL_PROFILE("amr-wind::ib::IB::post_init_actions");
    m_ib_levelset.setVal(1.0e30_rt);

    for (auto& ib : m_ibs) {
        ib->init_ib();
    }
}

void IB::post_regrid_actions() {}

void IB::pre_advance_work()
{
    BL_PROFILE("amr-wind::ib::IB::pre_advance_work");
}

void IB::pre_pressure_correction_work()
{
    BL_PROFILE("amr-wind::ib::IB::pre_pressure_correction_work");
    update_velocities();
}

void IB::post_pressure_correction_work()
{
    BL_PROFILE("amr-wind::ib::IB::pre_pressure_correction_work");
    update_velocities();
}

/** Update immersed boundary positions.
 *
 *  This method loops over all immersed boundaries and updates
 *  their position. Their new location may be prescribed or computed
 *
 *  \sa IB::update_velocities
 */
void IB::update_positions()
{
    BL_PROFILE("amr-wind::ib::IB::update_positions");
}

/** Provide updated velocities to immersed boundary instances
 *
 *  \sa IB::update_positions
 */
void IB::update_velocities()
{
    BL_PROFILE("amr-wind::ib::IB::update_velocity");
    for (auto& ib : m_ibs) {
        ib->update_velocities();
    }
}

/** Helper method to compute forces on all immersed boundary components
 */
void IB::compute_forces()
{
    BL_PROFILE("amr-wind::ib::IB::compute_forces");
    for (auto& ib : m_ibs) {
        ib->compute_forces();
    }
}

void IB::prepare_outputs()
{
    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string out_dir_prefix = post_dir + "/immersed_boundary";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }
    for (auto& ib : m_ibs) {
        ib->prepare_outputs(sname);
    }
}

void IB::post_advance_work()
{
    BL_PROFILE("amr-wind::ib::IB::post_advance_work");
    for (auto& ib : m_ibs) {
        ib->compute_forces();
        ib->update_positions();
        ib->write_outputs();
    }
}

} // namespace amr_wind::ib
