#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/overset/OversetManager.H"
#include "amr-wind/core/ExtSolver.H"
#include "amr-wind/wind_energy/ABL.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

CFDSim::CFDSim(amrex::AmrCore& mesh)
    : m_mesh(mesh)
    , m_repo(m_mesh)
    , m_pde_mgr(*this)
    , m_io_mgr(new IOManager(*this))
    , m_post_mgr(new PostProcessManager(*this))
    , m_ext_solver_mgr(new ExtSolverMgr)
    , m_helics(new HelicsStorage(*this))
{}

CFDSim::~CFDSim() = default;

void CFDSim::create_transport_model()
{
    std::string transport_model_name = "ConstTransport";
    {
        amrex::ParmParse pp("transport");
        pp.query("model", transport_model_name);
    }
    m_transport =
        transport::TransportModel::create(transport_model_name, *this);
}

void CFDSim::create_turbulence_model()
{
    std::string transport_model_name = "ConstTransport";
    std::string turbulence_model_name = "Laminar";
    {
        amrex::ParmParse pp("transport");
        pp.query("model", transport_model_name);
    }
    {
        amrex::ParmParse pp("turbulence");
        pp.query("model", turbulence_model_name);
    }

    const std::string identifier =
        turbulence_model_name + "-" + transport_model_name;
    m_turbulence = turbulence::TurbulenceModel::create(identifier, *this);
    m_turbulence->parse_model_coeffs();
}

void CFDSim::init_physics()
{
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> phys_names;
    pp.queryarr("physics", phys_names);

    for (const auto& phy : phys_names) {
        m_physics_mgr.create(phy, *this);
    }
}

void CFDSim::init_field_boundaries()
{
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> fb_names;
    pp.queryarr("field_boundaries", fb_names);

    for (const auto& fb : fb_names) {
        m_field_boundary_mgr.create(fb, *this);
    }
}

void CFDSim::activate_overset()
{
    amrex::ParmParse pp("overset");
    std::string otype = "TIOGA";

    m_overset_mgr = OversetManager::create(otype, *this);
}

bool CFDSim::has_overset() const { return (static_cast<bool>(m_overset_mgr)); }

void CFDSim::activate_mesh_map()
{
    amrex::ParmParse pp("geometry");
    std::string mesh_map_name; // default
    m_mesh_mapping = static_cast<bool>(pp.query("mesh_mapping", mesh_map_name));
    if (m_mesh_mapping) {
        m_mesh_map = MeshMap::create(mesh_map_name);
        m_mesh_map->declare_mapping_fields(*this, m_pde_mgr.num_ghost_state());
    }
}

bool CFDSim::is_anelastic() const
{
    if (m_physics_mgr.contains("ABL")) {
        const auto& abl = m_physics_mgr.template get<ABL>();
        return abl.anelastic().is_anelastic();
    }
    return false;
}
} // namespace amr_wind
