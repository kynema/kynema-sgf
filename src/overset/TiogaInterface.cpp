#include "amr-wind/overset/TiogaInterface.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/overset/overset_ops_routines.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_ParmParse.H"

#include <memory>
#include <numeric>
namespace amr_wind {

AMROversetInfo::AMROversetInfo(const int nglobal, const int nlocal)
    : level(nglobal)
    , mpi_rank(nglobal)
    , local_id(nglobal)
    , ilow(static_cast<size_t>(AMREX_SPACEDIM) * nglobal)
    , ihigh(static_cast<size_t>(AMREX_SPACEDIM) * nglobal)
    , dims(static_cast<size_t>(AMREX_SPACEDIM) * nglobal)
    , xlo(static_cast<size_t>(AMREX_SPACEDIM) * nglobal)
    , dx(static_cast<size_t>(AMREX_SPACEDIM) * nglobal)
    , global_idmap(nlocal)
    , iblank_node(nlocal)
    , iblank_cell(nlocal)
    , qcell(nlocal)
    , qnode(nlocal)
    , ngrids_global(nglobal)
    , ngrids_local(nlocal)
{}

// clang-format off
TiogaInterface::TiogaInterface(CFDSim& sim)
    : m_sim(sim)
    , m_iblank_cell(sim.repo().declare_int_field(
          "iblank_cell", 1, sim.pde_manager().num_ghost_state()))
    , m_iblank_node(sim.repo().declare_int_field(
          "iblank_node", 1, sim.pde_manager().num_ghost_state(), 1,
          FieldLoc::NODE))
    , m_mask_cell(sim.repo().declare_int_field(
          "mask_cell", 1, sim.pde_manager().num_ghost_state()))
    , m_mask_node(sim.repo().declare_int_field(
          "mask_node", 1, sim.pde_manager().num_ghost_state(), 1,
          FieldLoc::NODE))
{
    m_sim.io_manager().register_output_int_var(m_iblank_cell.name());
}

// clang-format on

void TiogaInterface::post_init_actions()
{
    amr_to_tioga_mesh();

    // Initialize masking so that all cells are active in solvers
    m_mask_cell.setVal(1);
    m_mask_node.setVal(1);
}

void TiogaInterface::post_regrid_actions()
{
    amr_to_tioga_mesh();
    amr_to_tioga_iblank();

    // Initialize masking so that all cells are active in solvers
    m_mask_cell.setVal(1);
    m_mask_node.setVal(1);
}

void TiogaInterface::pre_overset_conn_work()
{

    const auto& repo = m_sim.repo();
    const int num_ghost = m_sim.pde_manager().num_ghost_state();
    m_iblank_cell_host = repo.create_int_scratch_field_on_host(
        "iblank_cell_host", 1, num_ghost, FieldLoc::CELL);
    m_iblank_node_host = repo.create_int_scratch_field_on_host(
        "iblank_node_host", 1, num_ghost, FieldLoc::NODE);

    amr_to_tioga_iblank();

    m_iblank_cell.setVal(1);
    m_iblank_node.setVal(1);

    (*m_iblank_cell_host).setVal(1);
    (*m_iblank_node_host).setVal(1);
}

void TiogaInterface::post_overset_conn_work()
{

    const auto& repo = m_sim.repo();
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        htod_memcpy(m_iblank_cell(lev), (*m_iblank_cell_host)(lev), 0, 0, 1);
        htod_memcpy(m_iblank_node(lev), (*m_iblank_node_host)(lev), 0, 0, 1);

        m_iblank_cell(lev).FillBoundary(m_sim.mesh().Geom()[lev].periodicity());
        m_iblank_node(lev).FillBoundary(m_sim.mesh().Geom()[lev].periodicity());
    }

    overset_ops::iblank_to_mask(m_iblank_cell, m_mask_cell);
    overset_ops::iblank_to_mask(m_iblank_node, m_mask_node);

    // Update equation systems after a connectivity update
    m_sim.pde_manager().icns().post_regrid_actions();
    for (auto& eqn : m_sim.pde_manager().scalar_eqns()) {
        eqn->post_regrid_actions();
    }
    m_iblank_cell_host.reset();
    m_iblank_node_host.reset();
}

void TiogaInterface::register_solution(
    const std::vector<std::string>& cell_vars,
    const std::vector<std::string>& node_vars)
{
    auto& repo = m_sim.repo();
    const auto comp_counter =
        [&repo](int total, const std::string& fname) -> int {
        return total + repo.get_field(fname).num_comp();
    };
    const int ncell_vars =
        std::accumulate(cell_vars.begin(), cell_vars.end(), 0, comp_counter);
    const int nnode_vars =
        std::accumulate(node_vars.begin(), node_vars.end(), 0, comp_counter);
    const int num_ghost = m_sim.pde_manager().num_ghost_state();
    m_qcell = repo.create_scratch_field(ncell_vars, num_ghost, FieldLoc::CELL);
    m_qnode = repo.create_scratch_field(nnode_vars, num_ghost, FieldLoc::NODE);

    m_qcell_host = repo.create_scratch_field_on_host(
        ncell_vars, num_ghost, FieldLoc::CELL);
    m_qnode_host = repo.create_scratch_field_on_host(
        nnode_vars, num_ghost, FieldLoc::NODE);
    // Store field variable names for use in update_solution step
    m_cell_vars = cell_vars;
    m_node_vars = node_vars;

    const amrex::Real time_fillpatch = m_sim.is_during_overset_advance()
                                           ? m_sim.time().new_time()
                                           : m_sim.time().current_time();

    // Move cell variables into scratch field
    {
        int icomp = 0;
        for (const auto& cvar : m_cell_vars) {
            auto& fld = repo.get_field(cvar);
            const int ncomp = fld.num_comp();
            fld.fillpatch(time_fillpatch);
            field_ops::copy(*m_qcell, fld, 0, icomp, ncomp, num_ghost);
            icomp += ncomp;
        }
    }

    // Move cell variables into host scratch field
    {
        int icomp = 0;
        for (const auto& cvar : m_cell_vars) {
            auto& fld = repo.get_field(cvar);
            const int ncomp = fld.num_comp();
            fld.fillpatch(time_fillpatch);
            // Device to host copy happens here
            const int nlevels = repo.num_active_levels();
            for (int lev = 0; lev < nlevels; ++lev) {
                dtoh_memcpy((*m_qcell_host)(lev), fld(lev), 0, icomp, ncomp);
            }

            icomp += ncomp;
        }
    }
    // Move node variables into scratch field
    {
        int icomp = 0;
        for (const auto& cvar : m_node_vars) {
            auto& fld = repo.get_field(cvar);

            const int ncomp = fld.num_comp();
            fld.fillpatch(time_fillpatch);
            field_ops::copy(*m_qnode, fld, 0, icomp, ncomp, num_ghost);
            icomp += ncomp;
        }
        AMREX_ASSERT(nnode_vars == icomp);
    }
    // Copy node data from device to host scratch field
    {
        int icomp = 0;
        for (const auto& cvar : m_node_vars) {
            auto& fld = repo.get_field(cvar);
            const int ncomp = fld.num_comp();
            fld.fillpatch(time_fillpatch);
            // Device to host copy happens here
            const int nlevels = repo.num_active_levels();
            for (int lev = 0; lev < nlevels; ++lev) {
                dtoh_memcpy((*m_qnode_host)(lev), fld(lev), 0, icomp, ncomp);
            }
            icomp += ncomp;
        }
        AMREX_ASSERT(nnode_vars == icomp);
    }

    // Update data pointers for TIOGA exchange
    {
        int ilp = 0;
        const int nlevels = m_sim.repo().num_active_levels();
        auto& ad = *m_amr_data;
        amrex::Vector<amrex::Real*> qcellPtr(ad.qcell.size());
        amrex::Vector<amrex::Real*> qnodePtr(ad.qnode.size());
        for (int lev = 0; lev < nlevels; ++lev) {
            auto& qcfab = (*m_qcell)(lev);
            auto& qnfab = (*m_qnode)(lev);

            auto& qcfab_host = (*m_qcell_host)(lev);
            auto& qnfab_host = (*m_qnode_host)(lev);

            for (amrex::MFIter mfi(qcfab); mfi.isValid(); ++mfi) {
                // Copy from device to host
                ad.qcell.h_view[ilp] = qcfab_host[mfi].dataPtr();
                // Pass pointers to vector prior to d_view copy
                qcellPtr[ilp] = qcfab[mfi].dataPtr();
                // Copy from device to host
                ad.qnode.h_view[ilp] = qnfab_host[mfi].dataPtr();
                // Pass pointers to vector prior to d_view copy
                qnodePtr[ilp] = qnfab[mfi].dataPtr();

                ++ilp;
            }
        }
        // Host to device copy from standard vector to d_view
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, qcellPtr.begin(), qcellPtr.end(),
            ad.qcell.d_view.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, qnodePtr.begin(), qnodePtr.end(),
            ad.qnode.d_view.begin());
    }
}

void TiogaInterface::update_solution()
{
    auto& repo = m_sim.repo();

    const amrex::Real time_fillpatch = m_sim.is_during_overset_advance()
                                           ? m_sim.time().new_time()
                                           : m_sim.time().current_time();

    // Update cell variables on device
    {
        int icomp = 0;
        for (const auto& cvar : m_cell_vars) {
            auto& fld = repo.get_field(cvar);
            const int ncomp = fld.num_comp();
            // Host to device copy happens here
            const int nlevels = repo.num_active_levels();
            for (int lev = 0; lev < nlevels; ++lev) {
                htod_memcpy(fld(lev), (*m_qcell_host)(lev), icomp, 0, ncomp);
            }
            fld.fillpatch(time_fillpatch);
            icomp += ncomp;
        }
    }

    // Update nodal variables on device
    {
        int icomp = 0;
        for (const auto& cvar : m_node_vars) {
            auto& fld = repo.get_field(cvar);
            const int ncomp = fld.num_comp();
            // Host to device copy happens here
            const int nlevels = repo.num_active_levels();
            for (int lev = 0; lev < nlevels; ++lev) {
                htod_memcpy(fld(lev), (*m_qnode_host)(lev), icomp, 0, ncomp);
            }
            fld.fillpatch(time_fillpatch);
            icomp += ncomp;
        }
    }

    // Release memory to avoid holding onto scratch fields across regrids
    m_qcell.reset();
    m_qnode.reset();

    m_qcell_host.reset();
    m_qnode_host.reset();
}

void TiogaInterface::amr_to_tioga_mesh()
{
    BL_PROFILE("amr-wind::TiogaInterface::amr_to_tioga_mesh");
    auto& mesh = m_sim.mesh();
    const int nlevels = mesh.finestLevel() + 1;
    const int iproc = amrex::ParallelDescriptor::MyProc();
    const int nproc = amrex::ParallelDescriptor::NProcs();
    const auto* problo = mesh.Geom(0).ProbLo();

    int ngrids_global = 0;
    int ngrids_local = 0;
    for (int lev = 0; lev < nlevels; ++lev) {
        ngrids_global += static_cast<int>(mesh.boxArray(lev).size());

        const auto& dmap = mesh.DistributionMap(lev);
        AMREX_ALWAYS_ASSERT(
            dmap.size() <
            static_cast<amrex::Long>(std::numeric_limits<int>::max()));
        for (int d = 0; d < static_cast<int>(dmap.size()); ++d) {
            if (dmap[d] == iproc) {
                ++ngrids_local;
            }
        }
    }
    m_amr_data = std::make_unique<AMROversetInfo>(ngrids_global, ngrids_local);
    std::vector<int> lgrid_id(nproc, 0);

    int igp = 0; // Global index of the grid
    int ilp = 0; // Counter for local patches
    int iix = 0; // Index into the integer array

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ba = mesh.boxArray(lev);
        const auto& dm = mesh.DistributionMap(lev);
        const amrex::Real* dx = mesh.Geom(lev).CellSize();

        auto& ad = *m_amr_data;
        for (int d = 0; d < static_cast<int>(dm.size()); ++d) {
            ad.level.h_view[iix] = lev;      // AMR Level of this patch
            ad.mpi_rank.h_view[iix] = dm[d]; // MPI rank of this patch
            ad.local_id.h_view[iix] =
                lgrid_id[dm[d]]; // Local ID for this patch

            const auto& bx = ba[d];
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();

            const int ioff = AMREX_SPACEDIM * iix;
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                ad.ilow.h_view[ioff + i] = lo[i];
                ad.ihigh.h_view[ioff + i] = hi[i];
                ad.dims.h_view[ioff + i] = (hi[i] - lo[i]) + 1;

                ad.xlo.h_view[ioff + i] = problo[i] + (lo[i] * dx[i]);
                ad.dx.h_view[ioff + i] = dx[i];
            }

            if (iproc == dm[d]) {
                ad.global_idmap.h_view[ilp++] = igp;
            }

            // Increment array counter
            ++iix;
            // Increment global ID counter
            ++igp;
            // Increment local index
            ++lgrid_id[dm[d]];
        }
    }

    // Synchronize data on host/device
    m_amr_data->level.copy_to_device();
    m_amr_data->mpi_rank.copy_to_device();
    m_amr_data->local_id.copy_to_device();
    m_amr_data->ilow.copy_to_device();
    m_amr_data->ihigh.copy_to_device();
    m_amr_data->dims.copy_to_device();
    m_amr_data->xlo.copy_to_device();
    m_amr_data->dx.copy_to_device();
    m_amr_data->global_idmap.copy_to_device();
}
void TiogaInterface::amr_to_tioga_iblank()
{
    BL_PROFILE("amr-wind::TiogaInterface::amr_to_tioga_iblank");
    const auto& mesh = m_sim.mesh();
    const int nlevels = mesh.finestLevel() + 1;

    // Reset local patch counter
    int ilp = 0;
    auto& ibcell = m_sim.repo().get_int_field("iblank_cell");
    auto& ibnode = m_sim.repo().get_int_field("iblank_node");

    // Create standard vector of integer pointers here
    auto& ad = *m_amr_data;
    amrex::Vector<int*> tmpdataPtr(ad.iblank_cell.size());
    amrex::Vector<int*> tmpdataPtrn(ad.iblank_node.size());
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& ibfab = ibcell(lev);
        auto& ibnodefab = ibnode(lev);
        auto& ibfab_host = (*m_iblank_cell_host)(lev);
        auto& ibnodefab_host = (*m_iblank_node_host)(lev);

        for (amrex::MFIter mfi(ibfab); mfi.isValid(); ++mfi) {
            auto& ib = ibfab[mfi];
            auto& ibn = ibnodefab[mfi];
            auto& ib_host = ibfab_host[mfi];
            auto& ibn_host = ibnodefab_host[mfi];

            ad.iblank_cell.h_view[ilp] = ib_host.dataPtr();
            ad.iblank_node.h_view[ilp] = ibn_host.dataPtr();

            tmpdataPtr[ilp] = ib.dataPtr();
            tmpdataPtrn[ilp] = ibn.dataPtr();

            ++ilp;
        }
    }

    // Host to device copy from standard vector to
    // ad.iblank_cell.d_view
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tmpdataPtr.begin(), tmpdataPtr.end(),
        ad.iblank_cell.d_view.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tmpdataPtrn.begin(), tmpdataPtrn.end(),
        ad.iblank_node.d_view.begin());
}
} // namespace amr_wind
