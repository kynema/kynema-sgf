#ifndef FORESTDRAG_H
#define FORESTDRAG_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/integrals.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::forestdrag {

struct Forest
{
    int m_id{-1};
    amrex::Real m_type_forest;
    amrex::Real m_x_forest;
    amrex::Real m_y_forest;
    amrex::Real m_height_forest;
    amrex::Real m_diameter_forest;
    amrex::Real m_cd_forest;
    amrex::Real m_lai_forest;
    amrex::Real m_laimax_forest;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real lm() const
    {
        amrex::Real treelaimax = 0.0_rt;
        if (m_type_forest == 2) {
            const amrex::Real captured_height_forest = m_height_forest;
            const amrex::Real treeZm = m_laimax_forest * captured_height_forest;
            const int n = 100;
            const amrex::Real int_0_fh = amr_wind::utils::trapz(
                0, captured_height_forest - constants::TIGHT_TOL, n,
                [=](const amrex::Real x) {
                    const amrex::Real ratio =
                        (captured_height_forest - treeZm) /
                        (captured_height_forest - x);
                    const auto exponent = (x < treeZm) ? 6.0_rt : 0.5_rt;
                    return std::pow(ratio, exponent) *
                           std::exp(exponent * (1 - ratio));
                });
            treelaimax = m_lai_forest / int_0_fh;
        }
        return treelaimax;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    area_fraction(const amrex::Real z, const amrex::Real treelaimax) const
    {
        amrex::Real af = 0.0_rt;
        if (m_type_forest == 1) {
            af = m_lai_forest / m_height_forest;
        } else if (m_type_forest == 2) {
            const auto treeZm = m_laimax_forest * m_height_forest;
            const auto ratio =
                (m_height_forest - treeZm) / (m_height_forest - z);
            const auto exponent = (z < treeZm) ? 6.0_rt : 0.5_rt;
            af = treelaimax * std::pow(ratio, exponent) *
                 std::exp(exponent * (1 - ratio));
        }
        return af;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::RealBox real_bounding_box(
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo) const
    {
        const amrex::Real search_tol = 1.1_rt;
        const amrex::Real search_radius =
            0.5_rt * search_tol * m_diameter_forest;
        const auto x0 = m_x_forest - search_radius;
        const auto y0 = m_y_forest - search_radius;
        const auto z0 = prob_lo[2];
        const auto x1 = m_x_forest + search_radius;
        const auto y1 = m_y_forest + search_radius;
        const auto z1 = m_height_forest;
        return {x0, y0, z0, x1, y1, z1};
    }

    [[nodiscard]] amrex::Box bounding_box(const amrex::Geometry& geom) const
    {
        return utils::realbox_to_box(
            real_bounding_box(geom.ProbLoArray()), geom);
    }
};

/** Forestdrag Flow physics
 *  \ingroup physics
 */

class ForestDrag : public Physics::Register<ForestDrag>
{
public:
    static std::string identifier() { return "ForestDrag"; }

    explicit ForestDrag(CFDSim& sim);

    ~ForestDrag() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

    [[nodiscard]] amrex::Vector<Forest> read_forest(int level) const;

private:
    CFDSim& m_sim;
    Field& m_forest_drag;
    Field& m_forest_id;
    std::string m_forest_file{"forest.amrwind"};
};
} // namespace amr_wind::forestdrag

#endif
