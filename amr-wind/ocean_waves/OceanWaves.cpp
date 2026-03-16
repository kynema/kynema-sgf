#include "amr-wind/ocean_waves/OceanWaves.H"
#include "amr-wind/ocean_waves/OceanWavesModel.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/utilities/IOManager.H"

#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves {

OceanWaves::OceanWaves(CFDSim& sim)
    : m_sim(sim)
    , m_ow_levelset(sim.repo().declare_field("ow_levelset", 1, 3, 1))
    , m_ow_vof(sim.repo().declare_field("ow_vof", 1, 2, 1))
    , m_ow_velocity(
          sim.repo().declare_field("ow_velocity", AMREX_SPACEDIM, 3, 1))
{
    if (!sim.physics_manager().contains("MultiPhase")) {
        m_multiphase_mode = false;
    }
    m_ow_levelset.set_default_fillpatch_bc(sim.time());
    m_ow_vof.set_default_fillpatch_bc(sim.time());
    m_ow_velocity.set_default_fillpatch_bc(sim.time());

    m_ow_bndry = std::make_unique<OceanWavesBoundary>(sim);
}

OceanWaves::~OceanWaves() = default;

void OceanWaves::pre_init_actions()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::pre_init_actions");
    amrex::ParmParse pp(identifier());

    if (!(m_multiphase_mode ||
          m_sim.physics_manager().contains("TerrainDrag"))) {
        amrex::Abort(
            "OceanWaves requires MultiPhase or TerrainDrag physics to be "
            "active");
    }

    std::string label;
    pp.query("label", label);
    const std::string& tname = label;
    const std::string& prefix = identifier() + "." + tname;
    amrex::ParmParse pp1(prefix);

    std::string type;
    pp.query("type", type);
    pp1.query("type", type);
    AMREX_ALWAYS_ASSERT(!type.empty());

    m_owm = OceanWavesModel::create(type, m_sim, tname, 0);

    const std::string default_prefix = identifier() + "." + type;
    ::amr_wind::utils::MultiParser inp(default_prefix, prefix);

    m_owm->read_inputs(inp);
}

void OceanWaves::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::initialize_fields");
    m_owm->init_waves(level, geom, m_multiphase_mode);
}

void OceanWaves::post_init_actions()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::post_init_actions");
    m_ow_bndry->post_init_actions();
    m_owm->update_target_fields(m_sim.time().current_time());
    m_ow_bndry->record_boundary_data_time(m_sim.time().current_time());
    if (m_multiphase_mode) {
        m_owm->apply_relax_zones();
    }
    m_owm->reset_regrid_flag();
}

void OceanWaves::post_regrid_actions()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::post_regrid_actions");
    m_owm->record_regrid_flag();
}

void OceanWaves::pre_advance_work()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::pre_advance_work");
    // Update ow values for advection boundaries
    const amrex::Real adv_bdy_time =
        0.5_rt * (m_sim.time().current_time() + m_sim.time().new_time());
    m_owm->update_target_fields(adv_bdy_time);
    m_ow_bndry->record_boundary_data_time(adv_bdy_time);
}

void OceanWaves::pre_predictor_work()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::pre_predictor_work");
    // Update ow values for boundary fills at new time
    const amrex::Real bdy_fill_time = m_sim.time().new_time();
    m_owm->update_target_fields(bdy_fill_time);
    m_ow_bndry->record_boundary_data_time(bdy_fill_time);
}

void OceanWaves::post_advance_work()
{
    BL_PROFILE("amr-wind::ocean_waves::OceanWaves::post_init_actions");
    if (m_multiphase_mode) {
        m_owm->apply_relax_zones();
    }
    m_owm->reset_regrid_flag();
}

void OceanWaves::prepare_outputs()
{
    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string out_dir_prefix = post_dir + "/ocean_waves";
    const std::string sname =
        amrex::Concatenate(out_dir_prefix, m_sim.time().time_index());
    if (!amrex::UtilCreateDirectory(sname, 0755)) {
        amrex::CreateDirectoryFailed(sname);
    }

    m_owm->prepare_outputs(sname);
}

} // namespace amr_wind::ocean_waves
