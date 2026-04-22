#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "AMReX_REAL.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

using namespace amrex::literals;

namespace amr_wind {
class CFDSim;
class Field;
class IntField;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = FieldRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field value
 *  # tagging/t1.field_error = 10.0_rt 10.0_rt

 *  # Tagging based on local gradient
 *  tagging/t1.grad_error = 0.1_rt 0.1_rt 0.1_rt
 *  ```
 */
class FieldRefinement : public RefinementCriteria::Register<FieldRefinement>
{
public:
    static std::string identifier() { return "FieldRefinement"; }

    explicit FieldRefinement(const CFDSim& sim);

    ~FieldRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    template <typename MF>
    void tag(const int level, amrex::TagBoxArray& tags, const MF& mfab)
    {
        const bool tag_field = level <= m_max_lev_field;
        const bool tag_grad = level <= m_max_lev_grad;
        const auto& geom = m_sim.repo().mesh().Geom(level);
        const auto& prob_lo = geom.ProbLoArray();
        const auto& dx = geom.CellSizeArray();
        const auto tagging_box = m_tagging_box;

        const auto& tag_arrs = tags.arrays();
        const auto& farrs = mfab.const_arrays();
        if (tag_field) {
            const auto fld_err = m_field_error[level];
            amrex::ParallelFor(
                mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const amrex::RealVect coord = {AMREX_D_DECL(
                        prob_lo[0] + ((i + 0.5_rt) * dx[0]),
                        prob_lo[1] + ((j + 0.5_rt) * dx[1]),
                        prob_lo[2] + ((k + 0.5_rt) * dx[2]))};

                    if ((farrs[nbx](i, j, k) > fld_err) &&
                        (tagging_box.contains(coord))) {
                        tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
                    }
                });
        }

        if (tag_grad) {
            const auto gerr = m_grad_error[level];
            amrex::ParallelFor(
                mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    const amrex::RealVect coord = {AMREX_D_DECL(
                        prob_lo[0] + ((i + 0.5_rt) * dx[0]),
                        prob_lo[1] + ((j + 0.5_rt) * dx[1]),
                        prob_lo[2] + ((k + 0.5_rt) * dx[2]))};
                    const auto axp =
                        std::abs(farrs[nbx](i + 1, j, k) - farrs[nbx](i, j, k));
                    const auto ayp =
                        std::abs(farrs[nbx](i, j + 1, k) - farrs[nbx](i, j, k));
                    const auto azp =
                        std::abs(farrs[nbx](i, j, k + 1) - farrs[nbx](i, j, k));
                    const auto axm =
                        std::abs(farrs[nbx](i - 1, j, k) - farrs[nbx](i, j, k));
                    const auto aym =
                        std::abs(farrs[nbx](i, j - 1, k) - farrs[nbx](i, j, k));
                    const auto azm =
                        std::abs(farrs[nbx](i, j, k - 1) - farrs[nbx](i, j, k));
                    const auto ax = amrex::max<amrex::Real>(axp, axm);
                    const auto ay = amrex::max<amrex::Real>(ayp, aym);
                    const auto az = amrex::max<amrex::Real>(azp, azm);
                    if ((amrex::max<amrex::Real>(ax, ay, az) >= gerr) &&
                        (tagging_box.contains(coord))) {
                        tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
                    }
                });
        }
    }

private:
    const CFDSim& m_sim;

    Field* m_field{nullptr};
    IntField* m_int_field{nullptr};

    amrex::Vector<amrex::Real> m_field_error;
    amrex::Vector<amrex::Real> m_grad_error;

    int m_max_lev_field{-1};
    int m_max_lev_grad{-1};
    amrex::RealBox m_tagging_box;
};

} // namespace amr_wind

#endif /* FIELDREFINEMENT_H */
