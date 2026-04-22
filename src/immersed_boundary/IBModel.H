#ifndef IBMODEL_H
#define IBMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/core/MultiParser.H"

namespace amr_wind {

class CFDSim;

namespace ib {

/** Abstract representation of an immersed boundary in the flow field.
 *
 *  \ingroup immersed_boundary
 *
 * This class is intended to represent a single object (e.g., floating
 * platform, tower complex terrain) that exists in the computational domain. A
 * collection of immersed boundary objects (e.g., multiple cylinders) is modeled
 * using the amr_wind::IB physics class that holds a collection of instances of
 * IBModel and acts on them.
 */
class ImmersedBoundaryModel
    : public Factory<
          ImmersedBoundaryModel,
          CFDSim&,
          const std::string&,
          const int>
{
public:
    static std::string base_identifier() { return "ImmersedBoundaryModel"; }

    ~ImmersedBoundaryModel() override = default;

    [[nodiscard]] virtual const std::string& label() const = 0;

    [[nodiscard]] virtual int id() const = 0;

    virtual void read_inputs(const ::amr_wind::utils::MultiParser&) = 0;

    virtual void init_ib() = 0;

    virtual void update_positions() = 0;

    virtual void update_velocities() = 0;

    virtual void compute_forces() = 0;

    virtual void prepare_outputs(const std::string&) = 0;

    virtual void write_outputs() = 0;
};

/** Concrete implementation of the IBModel for different geometry types.
 *
 *  \ingroup immersed_boundary
 *
 * \tparam GeomTrait A geometry type trait that defines the behavior of the
 * different actions (initialization, velocity update, force computation)
 * that are executed by the IBModel class through the amr_wind::ib physics
 * driver class.
 *
 */
template <typename GeomTrait>
class IBModel : public ImmersedBoundaryModel::Register<IBModel<GeomTrait>>
{
private:
    //! Instance that holds all data related to a particular immersed boundary
    //! type
    typename GeomTrait::DataType m_data;

    //! Outputs processor
    ops::ProcessOutputsOp<GeomTrait> m_out_op;

public:
    static std::string identifier() { return GeomTrait::identifier(); }

    /**
     *  \param sim CFD simulation environment instance
     *  \param label Unique string identifier for this instance
     *  \param id Unique integer tag for this instance
     */
    IBModel(CFDSim& sim, const std::string& label, const int id)
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
        ops::ReadInputsOp<GeomTrait>()(m_data, pp);
        m_out_op.read_io_options(pp);
    }

    void update_positions() override {}

    void update_velocities() override { ops::UpdateVelOp<GeomTrait>()(m_data); }

    void compute_forces() override { ops::ComputeForceOp<GeomTrait>()(m_data); }

    void prepare_outputs(const std::string& out_dir) override
    {
        m_out_op.prepare_outputs(out_dir);
    }

    void write_outputs() override { m_out_op.write_outputs(); }

    void init_ib() override { ops::InitDataOp<GeomTrait>()(m_data); }
};

} // namespace ib
} // namespace amr_wind

#endif /* IBMODEL_H */
