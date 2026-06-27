#ifndef EBDRAG_H
#define EBDRAG_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"

namespace kynema_sgf::ebdrag {

/** EBdrag Flow physics
 *  \ingroup physics
 */

class EBDrag : public Physics::Register<EBDrag>
{
public:
    static std::string identifier() { return "EBDrag"; }

    explicit EBDrag(CFDSim& sim);

    ~EBDrag() override = default;

    void initialize_eb_fields(
        int level,
        const amrex::Geometry& geom,
        const amrex::MultiFab& vfrac_mf);

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override {}

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    //! Blanking field for terrain or buildings
    Field& m_eb_blank;
};
} // namespace kynema_sgf::ebdrag

#endif
