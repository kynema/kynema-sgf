#ifndef VOF_OPS_H
#define VOF_OPS_H

#include "src/equation_systems/vof/vof.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParReduce.H"

namespace kynema_sgf::pde {

/** Field registration operator for VOF equation system
 *  \ingroup vof
 */
template <typename Scheme>
struct FieldRegOp<VOF, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<VOF, Scheme>(time, repo);

        auto& levelset = repo.declare_cc_field("levelset", 1, 3, 1);
        auto& curvature = repo.declare_cc_field("interface_curvature", 1, 1, 1);

        levelset.set_default_fillpatch_bc(sim.time());
        curvature.set_default_fillpatch_bc(sim.time());

        // Register fields for output/restart
        auto& iomgr = sim.io_manager();
        iomgr.register_io_var(fields.field.name());

        return fields;
    }

    CFDSim& sim;
};

/** Right-hand side (RHS) evaluator for VOF transport equation
 */
template <typename Scheme>
struct ComputeRHSOp<VOF, Scheme>
{
    explicit ComputeRHSOp(PDEFields& fields_in) : fields(fields_in) {}

    void predictor_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real /*unused*/,
        bool /*unused*/)
    {}

    void corrector_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real /*unused*/,
        bool /*unused*/)
    {}

    void improve_explicit_diff(const amrex::Real /*unused*/) {}

    // data members
    PDEFields& fields;
};

/** Implementation of post-processing actions for the VOF equation
 *  \ingroup vof
 */

template <>
struct PostSolveOp<VOF>
{
    PostSolveOp(CFDSim& sim, PDEFields& fields) : m_sim(sim), m_fields(fields)
    {}

    /** Extrapolate vof laterally into terrain-blanked cells
     *
     *  Blanked cells that border unblanked cells in x or y take the average vof
     *  of those unblanked neighbors. Blanked cells surrounded by blanked cells
     *  take the average of their laterally neighboring cells that have already
     *  been assigned, propagating the values inward one layer per sweep.
     */
    void extrapolate_vof_into_terrain()
    {
        auto& repo = m_sim.repo();
        const auto& blanking = repo.get_int_field("terrain_blank");
        auto& vof = m_fields.field;
        const int nlevels = repo.num_active_levels();

        constexpr amrex::Real vof_tiny = 1.0e-8_rt;

        // Two components used to ping-pong between sweeps so that values
        // assigned during a sweep are not used as sources within that sweep.
        // A negative value flags a cell that has not been assigned yet.
        auto scratch = repo.create_scratch_field(2, 1);

        for (int lev = 0; lev < nlevels; ++lev) {
            auto& buf = (*scratch)(lev);
            const auto periodicity = m_sim.mesh().Geom(lev).periodicity();
            const auto ngrow = buf.nGrowVect();
            const auto domain = m_sim.mesh().Geom(lev).Domain();
            amrex::iMultiFab level_mask;
            if (lev < nlevels - 1) {
                level_mask = amrex::makeFineMask(
                    m_sim.mesh().boxArray(lev),
                    m_sim.mesh().DistributionMap(lev),
                    m_sim.mesh().boxArray(lev + 1),
                    m_sim.mesh().refRatio(lev), 1, 0);
            } else {
                level_mask.define(
                    m_sim.mesh().boxArray(lev),
                    m_sim.mesh().DistributionMap(lev), 1, 0,
                    amrex::MFInfo());
                level_mask.setVal(1);
            }
            buf.setVal(-1.0_rt);

            const auto blank_arrs = blanking(lev).const_arrays();
            const auto mask_arrs = level_mask.const_arrays();
            {
                const auto buf_arrs = buf.arrays();
                const auto vof_arrs = vof(lev).const_arrays();
                amrex::ParallelFor(
                    buf, amrex::IntVect(0),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                        if (mask_arrs[nbx](i, j, k) == 1 &&
                            domain.contains(i, j, k) &&
                            blank_arrs[nbx](i, j, k) == 0) {
                            buf_arrs[nbx](i, j, k, 0) = vof_arrs[nbx](i, j, k);
                        }
                    });
                amrex::Gpu::streamSynchronize();
            }
            buf.FillBoundary(0, 1, ngrow, periodicity);

            int src = 0;
            while (true) {
                const int dst = 1 - src;
                const auto buf_arrs = buf.arrays();
                auto nfilled = amrex::ParReduce(
                    amrex::TypeList<amrex::ReduceOpSum>{},
                    amrex::TypeList<amrex::Long>{}, buf, amrex::IntVect(0),
                    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
                        -> amrex::GpuTuple<amrex::Long> {
                        const auto& bf = buf_arrs[nbx];
                        const amrex::Real cur = bf(i, j, k, src);
                        if (cur >= 0.0_rt) {
                            bf(i, j, k, dst) = cur;
                            return {0};
                        }
                        amrex::Real sum = 0.0_rt;
                        int cnt = 0;
                        for (int ni = -1; ni <= 1; ++ni) {
                            for (int nj = -1; nj <= 1; ++nj) {
                                const int ii = i + ni;
                                const int jj = j + nj;
                                if (!domain.contains(ii, jj, k)) {
                                    continue;
                                }
                                const amrex::Real nval = bf(ii, jj, k, src);
                                if (nval >= 0.0_rt && !(ni == 0 && nj == 0)) {
                                    sum += nval;
                                    ++cnt;
                                }
                            }
                        }
                        if (cnt > 1) {
                            bf(i, j, k, dst) =
                                sum / static_cast<amrex::Real>(cnt);
                            // Snap to VOF bounds
                            if (bf(i, j, k, dst) < vof_tiny) {
                                bf(i, j, k, dst) = 0.0_rt;
                            } else if (bf(i, j, k, dst) > 1.0_rt - vof_tiny) {
                                bf(i, j, k, dst) = 1.0_rt;
                            }
                            return {1};
                        }
                        bf(i, j, k, dst) = -1.0_rt;
                        return {0};
                    });
                amrex::ParallelDescriptor::ReduceLongSum(nfilled);
                buf.FillBoundary(dst, 1, ngrow, periodicity);
                src = dst;
                // No cell was assigned during this sweep; nothing left to do
                if (nfilled == 0) {
                    break;
                }
            }

            const auto buf_arrs = buf.const_arrays();
            const auto vof_arrs = vof(lev).arrays();
            amrex::ParallelFor(
                vof(lev), amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const amrex::Real val = buf_arrs[nbx](i, j, k, src);
                    if (mask_arrs[nbx](i, j, k) == 1 &&
                        blank_arrs[nbx](i, j, k) == 1 &&
                        val >= 0.0_rt) {
                        vof_arrs[nbx](i, j, k) = val;
                    }
                });
            amrex::Gpu::streamSynchronize();
        }

        vof.fillpatch(m_sim.time().new_time());
    }

    void operator()(const amrex::Real /*unused*/)
    {
        m_fields.field.fillpatch(m_sim.time().new_time());

        if (m_sim.repo().int_field_exists("terrain_blank")) {
            extrapolate_vof_into_terrain();
        }

        auto& multiphase = m_sim.physics_manager().get<MultiPhase>();
        multiphase.set_density_via_vof();

        // Calculate n+1/2 density
        multiphase.set_nph_density();
        // Convert advalpha to advrho
        multiphase.calculate_advected_facedensity();
    }

    CFDSim& m_sim;
    PDEFields& m_fields;
};

} // namespace kynema_sgf::pde

#endif /* VOF_OPS_H */
