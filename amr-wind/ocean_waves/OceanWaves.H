#ifndef OCEAN_WAVES_H
#define OCEAN_WAVES_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/ocean_waves/boundary_ops/OceanWavesBoundary.H"

/** \defgroup ocean_waves Ocean Waves module
 *  Ocean Waves class
 */

namespace amr_wind {

class Field;

namespace ocean_waves {

class OceanWavesModel;

/** ocean waves modeling
 *
 *  \ingroup ocean_waves
 *
 *  This class provides an interface to model
 *
 *  \sa WaveModel
 */
class OceanWaves : public Physics::Register<OceanWaves>
{
public:
    static std::string identifier() { return "OceanWaves"; }

    explicit OceanWaves(CFDSim& /*sim*/);

    ~OceanWaves() override;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void pre_predictor_work() override;

    void post_advance_work() override;

    [[nodiscard]] const OceanWavesBoundary& ow_bndry() const
    {
        return *m_ow_bndry;
    }

protected:
    virtual void prepare_outputs();

private:
    CFDSim& m_sim;

    //! Unique pointer to the ocean waves model
    std::unique_ptr<OceanWavesModel> m_owm;

    //! Ocean waves boundary instance
    std::unique_ptr<OceanWavesBoundary> m_ow_bndry;

    //! Ocean waves target free surface levelset function
    Field& m_ow_levelset;

    //! Ocean waves target volume-of-fluid
    Field& m_ow_vof;

    //! Ocean waves target velocity
    Field& m_ow_velocity;

    //! Multiphase mode forces the velocity and vof in the domain
    //! When off, single-phase mode simply provides the target variables, and
    //! forcing is handled in other parts of the code
    bool m_multiphase_mode{true};
};

} // namespace ocean_waves
} // namespace amr_wind

#endif /* OCEAN_WAVES_H */
