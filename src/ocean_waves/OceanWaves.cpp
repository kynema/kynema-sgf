#include "src/ocean_waves/OceanWaves.H"
#include "src/ocean_waves/OceanWavesModel.H"
#include "src/CFDSim.H"
#include "src/core/FieldRepo.H"
#include "src/core/MultiParser.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "src/utilities/IOManager.H"

#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::ocean_waves {

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

    // Get the field boundaries named explicitly
    amrex::ParmParse pp_incflo("incflo");
    amrex::Vector<std::string> fb_names;
    pp_incflo.queryarr("field_boundaries", fb_names);
    bool bp_present = false;
    bool mpl_present = false;
    bool owb_present = false;
    for (const auto& fb : fb_names) {
        if (fb == "BoundaryPlane") {
            bp_present = true;
        }
        if (fb == "ModulatedPowerLaw") {
            mpl_present = true;
        }
        if (fb == "OceanWavesBoundary") {
            owb_present = true;
        }
    }

    // Add OceanWaves field boundary if not present unless conflicts are
    // detected
    int need_changes = static_cast<int>(!owb_present);
    if (bp_present) {
        // Check for input mode
        int pp_io_mode = -1;
        amrex::ParmParse pp_abl("ABL");
        pp_abl.query("bndry_io_mode", pp_io_mode);
        amrex::ParmParse pp_bdy("BoundaryPlane");
        pp_bdy.query("io_mode", pp_io_mode);
        if (pp_io_mode == 1) {
            // Turn off ow_bndry; will rely on bndry_plane for fills
            // Unless ow_bndry is not present, then no changes needed
            need_changes -= 1;
        }
    }
    if (mpl_present) {
        amrex::Abort(
            "OceanWavesBoundary: not currently compatible with Modulated Power "
            "Law implementation.");
    }

    if (need_changes == 1) {
        // Add OceanWavesBoundary to field boundaries list
        fb_names.push_back("OceanWavesBoundary");
    }
    if (need_changes == -1) {
        // Remove OceanWavesBoundary from field boundaries list and allow
        // BoundaryPlane to handle fills
        amrex::Print() << "OceanWavesBoundary: detected conflict with "
                       << "BoundaryPlane; removing OceanWavesBoundary from "
                       << "field boundaries list\n";
        fb_names.erase(
            std::remove(fb_names.begin(), fb_names.end(), "OceanWavesBoundary"),
            fb_names.end());
    }

    if (need_changes != 0) {
        // Update the input database with the new field boundaries list
        pp_incflo.addarr("field_boundaries", fb_names);
    }
}

OceanWaves::~OceanWaves() = default;

void OceanWaves::pre_init_actions()
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::pre_init_actions");
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
    ::kynema_sgf::utils::MultiParser inp(default_prefix, prefix);

    m_owm->read_inputs(inp);
}

void OceanWaves::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::initialize_fields");
    m_owm->init_waves(level, geom, m_multiphase_mode);
}

void OceanWaves::post_init_actions()
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::post_init_actions");
    m_owm->update_target_fields(m_sim.time().current_time());
    if (m_multiphase_mode) {
        m_owm->apply_relax_zones();
    }
    m_owm->reset_regrid_flag();
}

void OceanWaves::post_regrid_actions()
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::post_regrid_actions");
    m_owm->record_regrid_flag();
}

void OceanWaves::pre_advance_work()
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::pre_advance_work");
    // Update ow values for advection boundaries
    const amrex::Real adv_bdy_time =
        0.5_rt * (m_sim.time().current_time() + m_sim.time().new_time());
    m_owm->update_target_fields(adv_bdy_time);
}

void OceanWaves::pre_predictor_work()
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::pre_predictor_work");
    // Update ow values for boundary fills at new time
    const amrex::Real bdy_fill_time = m_sim.time().new_time();
    m_owm->update_target_fields(bdy_fill_time);
}

void OceanWaves::post_advance_work()
{
    BL_PROFILE("kynema-sgf::ocean_waves::OceanWaves::post_advance_work");
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

} // namespace kynema_sgf::ocean_waves
