#include <vector>

#include "src/overset/SGFTiogaIface.H"
#include "src/CFDSim.H"
#include "src/overset/TiogaInterface.H"
#include "TiogaMeshInfo.h"
#include "tioga.h"

namespace kynema_sgf {

namespace {

template <typename T1, typename T2>
void amr_to_tioga(T1& lhs, T2& rhs)
{
    lhs.sz = rhs.size();
    lhs.hptr = rhs.h_view.data();
    lhs.dptr = rhs.d_view.data();
}

} // namespace

SGFTiogaIface::SGFTiogaIface(kynema_sgf::CFDSim& sim, TIOGA::tioga& tg)
    : m_sim(sim), m_tg(tg), m_info(new TIOGA::AMRMeshInfo)
{}

SGFTiogaIface::~SGFTiogaIface() = default;

void SGFTiogaIface::pre_overset_conn_work()
{
    m_sim.overset_manager()->pre_overset_conn_work();
    register_mesh();
}

void SGFTiogaIface::post_overset_conn_work()
{
    m_sim.overset_manager()->post_overset_conn_work();
}

void SGFTiogaIface::register_mesh()
{
    BL_PROFILE("kynema_sgf::SGFTiogaIface::register_mesh");
    const int num_ghost = m_sim.pde_manager().num_ghost_state();

    auto* amr_tg_iface =
        dynamic_cast<kynema_sgf::TiogaInterface*>(m_sim.overset_manager());
    if (amr_tg_iface == nullptr) {
        amrex::Abort("Dynamic cast to TiogaInterface failed");
        return;
    }
    auto& ad = amr_tg_iface->amr_overset_info();
    auto& mi = *m_info;

    mi.ngrids_global = ad.ngrids_global;
    mi.ngrids_local = ad.ngrids_local;
    mi.num_ghost = num_ghost;
    amr_to_tioga(mi.level, ad.level);
    amr_to_tioga(mi.mpi_rank, ad.mpi_rank);
    amr_to_tioga(mi.local_id, ad.local_id);
    amr_to_tioga(mi.ilow, ad.ilow);
    amr_to_tioga(mi.ihigh, ad.ihigh);
    amr_to_tioga(mi.dims, ad.dims);
    amr_to_tioga(mi.xlo, ad.xlo);
    amr_to_tioga(mi.dx, ad.dx);
    amr_to_tioga(mi.global_idmap, ad.global_idmap);
    amr_to_tioga(mi.iblank_node, ad.iblank_node);
    amr_to_tioga(mi.iblank_cell, ad.iblank_cell);
    amr_to_tioga(mi.qcell, ad.qcell);
    amr_to_tioga(mi.qnode, ad.qnode);

    // Will reset when we register solution
    mi.nvar_cell = 0;
    mi.nvar_node = 0;

    m_tg.register_amr_grid(m_info.get());
}

void SGFTiogaIface::register_solution(
    const std::vector<std::string>& cell_vars,
    const std::vector<std::string>& node_vars)
{
    auto* amr_tg_iface =
        dynamic_cast<kynema_sgf::TiogaInterface*>(m_sim.overset_manager());
    if (amr_tg_iface == nullptr) {
        amrex::Abort("Dynamic cast to TiogaInterface failed");
        return;
    }
    amr_tg_iface->register_solution(cell_vars, node_vars);
    auto& qcell = amr_tg_iface->qvars_cell();
    auto& qnode = amr_tg_iface->qvars_node();

    // Ensure that ghost cells are consistent
    AMREX_ALWAYS_ASSERT(qcell.num_grow()[0] == qnode.num_grow()[0]);

    auto& ad = amr_tg_iface->amr_overset_info();
    auto& mi = *m_info;
    mi.nvar_cell = qcell.num_comp();
    mi.nvar_node = qnode.num_comp();
    amr_to_tioga(mi.qcell, ad.qcell);
    amr_to_tioga(mi.qnode, ad.qnode);
    m_tg.register_amr_solution();
}

void SGFTiogaIface::update_solution()
{
    m_sim.overset_manager()->update_solution();
}

} // namespace kynema_sgf
