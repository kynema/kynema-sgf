#ifndef WAVE_MODEL_H
#define WAVE_MODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/ocean_waves/OceanWavesTypes.H"
#include "amr-wind/ocean_waves/OceanWavesOps.H"
#include "amr-wind/core/MultiParser.H"

namespace amr_wind {

class CFDSim;

namespace ocean_waves {

/** Abstract representation of a wave model used for wave generation
 *
 *  \ingroup ocean_waves
 *
 * This class is intended to represent a single wave model
 */
class OceanWavesModel
    : public Factory<OceanWavesModel, CFDSim&, const std::string&, const int>
{
public:
    static std::string base_identifier() { return "OceanWaveModel"; }

    ~OceanWavesModel() override = default;

    [[nodiscard]] virtual const std::string& label() const = 0;

    [[nodiscard]] virtual int id() const = 0;

    virtual void read_inputs(const ::amr_wind::utils::MultiParser&) = 0;

    virtual void init_waves(int, const amrex::Geometry&, bool) = 0;

    virtual void update_target_fields(amrex::Real) = 0;

    virtual void apply_relax_zones() = 0;

    virtual void prepare_outputs(const std::string&) = 0;

    virtual void write_outputs() = 0;

    virtual void record_regrid_flag() = 0;

    virtual void reset_regrid_flag() = 0;
};

/** Concrete implementation of the WaveModel for different wave theories.
 *
 *  \ingroup ocean_waves
 *
 * \tparam WaveTheoryTrait A wave theory trait that defines the behavior of the
 * different actions (initialization, velocity update, free surface update)
 * that are executed by the WaveModel class through the amr_wind::ocean_waves
 * physics driver class.
 *
 */
template <typename WaveTheoryTrait>
class OWModel : public OceanWavesModel::Register<OWModel<WaveTheoryTrait>>
{
private:
    //! Instance that holds all data related to a particular wave type
    typename WaveTheoryTrait::DataType m_data;

    //! Outputs processor
    ops::ProcessOutputsOp<WaveTheoryTrait> m_out_op;

public:
    static std::string identifier() { return WaveTheoryTrait::identifier(); }

    /**
     *  \param sim CFD simulation environment instance
     *  \param label Unique string identifier for this instance
     *  \param id Unique integer tag for this instance
     */
    OWModel(CFDSim& sim, const std::string& label, const int id)
        : m_data(sim, label, id), m_out_op(m_data)
    {}

    //! Return the unique identifier (name) for this instance
    [[nodiscard]] const std::string& label() const override
    {
        return m_data.info().label;
    }

    //! Return the unique tag (integer ID) for this instance
    [[nodiscard]] int id() const override { return m_data.info().id; }

    void read_inputs(const ::amr_wind::utils::MultiParser& pp) override
    {
        ops::ReadInputsOp<WaveTheoryTrait>()(m_data, pp);
        m_out_op.read_io_options(pp);
    }

    void update_target_fields(const amrex::Real time) override
    {
        ops::UpdateTargetFieldsOp<WaveTheoryTrait>()(m_data, time);
        ops::UpdateTargetFieldsSecondaryOp<WaveTheoryTrait>()(m_data);
    }

    void apply_relax_zones() override
    {
        ops::ApplyRelaxZonesOp<WaveTheoryTrait>()(m_data);
    }

    void prepare_outputs(const std::string& out_dir) override
    {
        m_out_op.prepare_outputs(out_dir);
    }

    void record_regrid_flag() override { m_data.meta().regrid_occurred = true; }

    void reset_regrid_flag() override { m_data.meta().regrid_occurred = false; }

    void write_outputs() override { m_out_op.write_outputs(); }

    void init_waves(
        int level, const amrex::Geometry& geom, bool multiphase_mode) override
    {
        ops::InitDataOp<WaveTheoryTrait>()(
            m_data, level, geom, multiphase_mode);
    }
};

} // namespace ocean_waves
} // namespace amr_wind

#endif /* OWMODEL_H */
