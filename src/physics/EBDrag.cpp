#include "src/physics/EBDrag.H"
#include "src/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "src/utilities/trig_ops.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>

using namespace amrex::literals;

namespace kynema_sgf::ebdrag {

namespace {} // namespace

EBDrag::EBDrag(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_eb_blank(sim.repo().declare_field("eb_blank", 1, 1, 1))
    , m_eb_drag(sim.repo().declare_field("eb_drag", 1, 1, 1))
{
    amrex::ParmParse pp("ebdrag");
    m_sim.io_manager().register_output_var("eb_drag");
    m_sim.io_manager().register_output_var("eb_blank");
    m_eb_blank.setVal(0.0);
    m_eb_drag.setVal(0.0);
    amrex::ParmParse pp_eb("eb2");
    std::string geom;
    pp_eb.query("geom_type", geom);
    amrex::Print() << "[EB] Building geometry with geom_type = " << geom
                   << "\n";
}

void EBDrag::initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/)
{}

void EBDrag::initialize_eb_fields(
    int level, const amrex::Geometry& geom, const amrex::MultiFab& vfrac_mf)
{

    BL_PROFILE("kynema-sgf::" + this->identifier() + "::initialize_fields");
    // const auto& dx = geom.CellSizeArray();
    // const auto& prob_lo = geom.ProbLoArray();
    // const auto& prob_hi = geom.ProbHiArray();
    auto& blanking = m_eb_blank(level);
    auto& drag = m_eb_drag(level);
    const auto vfrac = vfrac_mf.arrays();
    auto levelBlanking = blanking.arrays();
    // auto levelDrag = drag.arrays();

    amrex::ParallelFor(
        blanking, m_eb_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            levelBlanking[nbx](i, j, k, 0) = vfrac[nbx](i, j, k);
        });
    amrex::Gpu::streamSynchronize();
}

void EBDrag::post_init_actions()
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::post_init_actions");
    const int nlevels = m_sim.repo().num_active_levels();
    int nghost = 1;
    amrex::Print() << "[EB] post_regrid_actions: rebuilding EB geometry and "
                      "updating blanking field on all levels\n";
    for (int lev = 0; lev < nlevels; ++lev) {
        const amrex::Geometry& geom = m_sim.repo().mesh().Geom(lev);
        const amrex::BoxArray& ba = m_sim.repo().mesh().boxArray(lev);
        const amrex::DistributionMapping& dm =
            m_sim.repo().mesh().DistributionMap(lev);
        amrex::EB2::Build(geom, lev, 0, 10, true);
        std::unique_ptr<amrex::EBFArrayBoxFactory> ebfactory = makeEBFabFactory(
            geom, ba, dm, {nghost, nghost, nghost}, amrex::EBSupport::full);
        const auto& vfrac_mf = ebfactory->getVolFrac();
        initialize_eb_fields(lev, m_sim.repo().mesh().Geom(lev), vfrac_mf);
        amrex::Print() << "[EB] EB2::Build complete. IndexSpace empty? "
                       << (amrex::EB2::IndexSpace::empty() ? "YES (problem!)"
                                                           : "NO (ok)")
                       << "\n";
    }
}

void EBDrag::pre_advance_work()
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::pre_advance_work");
}

void EBDrag::post_regrid_actions() {}

} // namespace kynema_sgf::ebdrag
