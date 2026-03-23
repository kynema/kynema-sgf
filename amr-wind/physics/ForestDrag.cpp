#include "amr-wind/physics/ForestDrag.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::forestdrag {

ForestDrag::ForestDrag(CFDSim& sim)
    : m_sim(sim)
    , m_forest_drag(sim.repo().declare_field("forest_drag", 1, 1, 1))
    , m_forest_id(sim.repo().declare_field("forest_id", 1, 1, 1))
{

    amrex::ParmParse pp(identifier());
    pp.query("forest_file", m_forest_file);

    m_sim.io_manager().register_output_var("forest_drag");
    m_sim.io_manager().register_output_var("forest_id");

    m_forest_drag.setVal(0.0_rt);
    m_forest_id.setVal(-1.0_rt);
    m_forest_id.set_default_fillpatch_bc(m_sim.time());
    m_forest_drag.set_default_fillpatch_bc(m_sim.time());
}

void ForestDrag::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    const auto forests = read_forest(level);
    amrex::Gpu::DeviceVector<Forest> d_forests(forests.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, forests.begin(), forests.end(),
        d_forests.begin());

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& drag = m_forest_drag(level);
    auto& fst_id = m_forest_id(level);
    drag.setVal(0.0_rt);
    fst_id.setVal(-1.0_rt);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(m_forest_drag(level)); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.growntilebox();
        for (int nf = 0; nf < static_cast<int>(forests.size()); nf++) {
            const auto bxi = vbx & forests[nf].bounding_box(geom);
            if (!bxi.isEmpty()) {
                const auto& levelDrag = drag.array(mfi);
                const auto& levelId = fst_id.array(mfi);
                const auto* forests_ptr = d_forests.data();
                amrex::ParallelFor(
                    bxi, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        const auto x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
                        const auto y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
                        const auto z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
                        const auto& fst = forests_ptr[nf];
                        const auto radius = std::sqrt(
                            ((x - fst.m_x_forest) * (x - fst.m_x_forest)) +
                            ((y - fst.m_y_forest) * (y - fst.m_y_forest)));
                        if (z <= fst.m_height_forest &&
                            radius <= (0.5_rt * fst.m_diameter_forest)) {
                            const auto treelaimax = fst.lm();
                            levelId(i, j, k) = fst.m_id;
                            levelDrag(i, j, k) +=
                                fst.m_cd_forest *
                                fst.area_fraction(z, treelaimax);
                        }
                    });
            }
        }
    }
}

void ForestDrag::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

amrex::Vector<Forest> ForestDrag::read_forest(const int level) const
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::read_forest");

    std::ifstream file(m_forest_file, std::ios::in);
    if (!file.good()) {
        amrex::Abort("Cannot find file " + m_forest_file);
    }

    //! TreeType xc yc height diameter cd lai laimax
    amrex::Vector<Forest> forests;
    const auto& geom = m_sim.repo().mesh().Geom(level);
    const auto& ba = m_sim.repo().mesh().boxArray(level);
    int cnt = 0;
    amrex::Real value1, value2, value3, value4, value5, value6, value7, value8;
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
        Forest f;
        f.m_id = cnt;
        f.m_type_forest = value1;
        f.m_x_forest = value2;
        f.m_y_forest = value3;
        f.m_height_forest = value4;
        f.m_diameter_forest = value5;
        f.m_cd_forest = value6;
        f.m_lai_forest = value7;
        f.m_laimax_forest = value8;

        // Only keep a forest if the rank owns it
        const auto bx = f.bounding_box(geom);
        if (ba.intersects(bx)) {
            forests.push_back(f);
        }
        cnt++;
    }
    file.close();
    return forests;
}
} // namespace amr_wind::forestdrag
