#include "amr-wind/physics/mms/MMS.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "masa.h"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::mms {

MMS::MMS(const CFDSim& sim)
    : m_time(sim.time())
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_mms_vel_source(sim.repo().declare_cc_field(
          "mms_vel_source", m_velocity.num_comp(), m_velocity.num_grow()[0], 1))
{
    std::string masa_name;
    amrex::ParmParse pp_mms("MMS");
    pp_mms.get("masa_name", masa_name);
    masa_init("mms", masa_name.c_str());
    for (auto& it : m_params_map) {
        pp_mms.query(it.first.c_str(), it.second);
    }
    amrex::ParmParse pp_transport("transport");
    pp_transport.query("viscosity", m_params_map["nu"]);
    for (const auto& it : m_params_map) {
        masa_set_param(it.first.c_str(), it.second);
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        masa_display_param();
        std::ofstream f;
        f.open(m_output_fname.c_str());
        f << std::setw(m_w) << "time" << std::setw(m_w) << "L2_u"
          << std::setw(m_w) << "L2_v" << std::setw(m_w) << "L2_w" << std::endl;
        f.close();
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::MMSFieldInit
 */
void MMS::initialize_fields(int level, const amrex::Geometry& geom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    amrex::MultiFab h_density(
        density.boxArray(), density.distributionMap, density.nComp(),
        density.nGrow(), amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab h_velocity(
        velocity.boxArray(), velocity.distributionMap, velocity.nComp(),
        velocity.nGrow(), amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(h_density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.tilebox();
        const auto& vel = h_velocity.array(mfi);
        const auto& den = h_density.array(mfi);

        amrex::LoopOnCpu(vbx, [=](int i, int j, int k) {
            const amrex::Real x = problo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = problo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = problo[2] + (k + 0.5_rt) * dx[2];

            den(i, j, k) = 1.0_rt;
            vel(i, j, k, 0) = masa_eval_3d_exact_u(x, y, z);
            vel(i, j, k, 1) = masa_eval_3d_exact_v(x, y, z);
            vel(i, j, k, 2) = masa_eval_3d_exact_w(x, y, z);
        });
    }
    amrex::htod_memcpy(density, h_density);
    amrex::htod_memcpy(velocity, h_velocity);
}

/** Fill the MMS source term.
 */
void MMS::fill_src()
{
    m_mms_vel_source.setVal(0.0_rt);

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& problo = m_mesh.Geom(lev).ProbLoArray();
        auto& mms_src_term = m_mms_vel_source(lev);
        amrex::MultiFab h_mms_src_term(
            mms_src_term.boxArray(), mms_src_term.distributionMap,
            mms_src_term.nComp(), mms_src_term.nGrow(),
            amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mms_src_term, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& mms_src = h_mms_src_term.array(mfi);

            amrex::LoopOnCpu(bx, [=](int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + 0.5_rt) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5_rt) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5_rt) * dx[2];

                mms_src(i, j, k, 0) = masa_eval_3d_source_u(x, y, z);
                mms_src(i, j, k, 1) = masa_eval_3d_source_v(x, y, z);
                mms_src(i, j, k, 2) = masa_eval_3d_source_w(x, y, z);
            });
        }
        amrex::htod_memcpy(mms_src_term, h_mms_src_term);
    }
}

void MMS::post_init_actions() { fill_src(); }

void MMS::post_regrid_actions() { fill_src(); }

amrex::Real
MMS::compute_error(const int comp, const Field& field, amr_wind::mms::FuncDef f)
{
    amrex::Real error = 0.0_rt;

    const int nlevels = m_repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev),
                m_mesh.boxArray(lev + 1), m_mesh.refRatio(lev), 1, 0);
        } else {
            level_mask.define(
                m_mesh.boxArray(lev), m_mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }
        amrex::iMultiFab h_level_mask(
            level_mask.boxArray(), level_mask.distributionMap,
            level_mask.nComp(), level_mask.nGrow(),
            amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        amrex::dtoh_memcpy(h_level_mask, level_mask);

        const auto& dx = m_mesh.Geom(lev).CellSizeArray();
        const auto& problo = m_mesh.Geom(lev).ProbLoArray();
        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        const auto& fld = field(lev);
        amrex::MultiFab h_fld(
            fld.boxArray(), fld.distributionMap, fld.nComp(), fld.nGrow(),
            amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        amrex::dtoh_memcpy(h_fld, fld);
        for (amrex::MFIter mfi(fld); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& field_arr = h_fld.array(mfi);
            const auto& mask_arr = h_level_mask.array(mfi);

            amrex::Real err_fab = 0.0_rt;
            amrex::LoopOnCpu(vbx, [=, &err_fab](int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + 0.5_rt) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5_rt) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5_rt) * dx[2];
                const amrex::Real u = field_arr(i, j, k, comp);
                const amrex::Real u_exact = f(x, y, z);
                err_fab += cell_vol * mask_arr(i, j, k) * (u - u_exact) *
                           (u - u_exact);
            });
            error += err_fab;
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(error);

    const amrex::Real total_vol = m_mesh.Geom(0).ProbDomain().volume();
    return std::sqrt(error / total_vol);
}

void MMS::post_advance_work()
{

    const amrex::Real u_mms_err =
        compute_error(0, m_velocity, masa_eval_3d_exact_u);
    const amrex::Real v_mms_err =
        compute_error(1, m_velocity, masa_eval_3d_exact_v);
    const amrex::Real w_mms_err =
        compute_error(2, m_velocity, masa_eval_3d_exact_w);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f;
        f.open(m_output_fname.c_str(), std::ios_base::app);
        f << std::setprecision(12) << std::setw(m_w) << m_time.current_time()
          << std::setw(m_w) << u_mms_err << std::setw(m_w) << v_mms_err
          << std::setw(m_w) << w_mms_err << std::endl;
        f.close();
    }
}

} // namespace amr_wind::mms
