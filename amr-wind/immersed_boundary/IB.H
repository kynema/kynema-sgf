#ifndef IB_H
#define IB_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/FieldRepo.H"

/** \defgroup immersed_boundary Immersed boundary module
 *  Immersed boundary class
 */

namespace amr_wind {

class Field;

namespace ib {

class ImmersedBoundaryModel;

/** Immersed boundary modeling for non-blade wind components and complex terrain
 *
 *  \ingroup immersed_boundary
 *
 *  This class provides an interface to model
 *
 *  \sa ImmersedBoundaryModel
 */
class IB : public Physics::Register<IB>
{
public:
    static std::string identifier() { return "IB"; }

    explicit IB(CFDSim& /*sim*/);

    ~IB() override;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

    void pre_pressure_correction_work() override;

    void post_pressure_correction_work() override;

protected:
    //! Total number of immersed boundaries

    [[nodiscard]] int num_ibs() const { return static_cast<int>(m_ibs.size()); }

    virtual void prepare_outputs();

private:
    static void update_positions();

    void update_velocities();

    void compute_forces();

    CFDSim& m_sim;

    std::vector<std::unique_ptr<ImmersedBoundaryModel>> m_ibs;

    //! Immersed boundary levelset for cell-centers
    Field& m_ib_levelset;

    //! Immersed boundary normal vector defined on cell-centers
    Field& m_ib_normal;
};

} // namespace ib
} // namespace amr_wind

#endif /* IB_H */
