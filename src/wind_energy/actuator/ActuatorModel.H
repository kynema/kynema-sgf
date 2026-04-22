#ifndef ACTUATORMODEL_H
#define ACTUATORMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/ActSrcLineOp.H"
#include "amr-wind/wind_energy/actuator/turbine/ActSrcDiskOp_Turbine.H"
#include "amr-wind/wind_energy/actuator/disk/ActSrcDiskOp.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind {

class CFDSim;

namespace actuator {

/** Abstract representation of an actuator source in the flow field.
 *
 *  \ingroup actuator
 *
 *  This class is intended to represent a single object (e.g., turbine, fixed
 *  wing) that exists in the computational domain. A collection of source
 *  objects (e.g., turbines in a wind farm) is modeled using the
 *  amr_wind::Actuator physics class that holds a collection of instances of
 *  ActuatorModel and acts on them.
 */
class ActuatorModel
    : public Factory<ActuatorModel, CFDSim&, const std::string&, const int>
{
public:
    static std::string base_identifier() { return "ActuatorModel"; }

    ~ActuatorModel() override = default;

    [[nodiscard]] virtual const std::string& label() const = 0;

    [[nodiscard]] virtual int id() const = 0;

    [[nodiscard]] virtual const ActInfo& info() const = 0;

    virtual ActInfo& info() = 0;

    virtual void read_inputs(const utils::ActParser&) = 0;

    virtual void determine_influenced_procs() = 0;

    virtual void determine_root_proc(amrex::Vector<int>&) = 0;

    virtual void init_actuator_source() = 0;

    [[nodiscard]] virtual int num_velocity_points() const = 0;

    virtual void update_positions(VecSlice&) = 0;

    virtual void update_fields(const VecSlice&, const RealSlice&) = 0;

    virtual void compute_forces() = 0;

    virtual void compute_source_term(
        int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom) = 0;

    virtual void prepare_outputs(const std::string&) = 0;

    virtual void write_outputs() = 0;
};

/** Concrete implementation of the ActuatorModel for different actuator types.
 *
 *  \ingroup actuator
 *
 *  \tparam ActTrait An actuator type trait that defines the behavior of the
 *  different actions (initialization, velocity sampling, force computation, and
 *  momentum source term computations) that are executed by the ActuatorModel
 *  class through the amr_wind::Actuator physics driver class.
 *
 *  \tparam SrcTrait A source type trait that defines how the forcing term is
 *  computed, e.g., if the turbine rotor is represented as a line or a disk. For
 *  most cases, like the fixed wing case, this will just be line, therefore, it
 *  is set as default for convenience.
 */
template <typename ActTrait, typename SrcTrait = ActSrcLine>
class ActModel : public ActuatorModel::Register<ActModel<ActTrait, SrcTrait>>
{
private:
    //! Instance that holds all data related to a particular actuator type
    typename ActTrait::DataType m_data;

    //! Source computation operator
    ops::ActSrcOp<ActTrait, SrcTrait> m_src_op;

    //! Outputs processor
    ops::ProcessOutputsOp<ActTrait, SrcTrait> m_out_op;

public:
    static std::string identifier()
    {
        return ActTrait::identifier() + SrcTrait::identifier();
    }

    /**
     *  \param sim CFD simulation environment instance
     *  \param label Unique string identifier for this instance
     *  \param id Unique integer tag for this instance
     */
    ActModel(CFDSim& sim, const std::string& label, const int id)
        : m_data(sim, label, id), m_src_op(m_data), m_out_op(m_data)
    {}

    //! Return the unique identifier (name) for this instance
    [[nodiscard]] const std::string& label() const override
    {
        return m_data.info().label;
    }

    //! Return the unique tag (integer ID) for this instance
    [[nodiscard]] int id() const override { return m_data.info().id; }

    //! Return the core info object for this actuator instance
    [[nodiscard]] const ActInfo& info() const override { return m_data.info(); }

    //! Return the core info object for this actuator instance
    ActInfo& info() override { return m_data.info(); }

    //! Return the meta info object for this actuator instance
    [[nodiscard]] const auto& meta() const { return m_data.meta(); }

    void read_inputs(const utils::ActParser& pp) override
    {
        ops::ReadInputsOp<ActTrait, SrcTrait>()(m_data, pp);
        m_out_op.read_io_options(pp);
    }

    void determine_influenced_procs() override
    {
        ops::determine_influenced_procs<ActTrait>(m_data);
    }

    void determine_root_proc(amrex::Vector<int>& act_proc_count) override;

    [[nodiscard]] int num_velocity_points() const override;

    void update_positions(VecSlice& pos) override
    {
        ops::UpdatePosOp<ActTrait, SrcTrait>()(m_data);
        const auto& vpos = m_data.grid().vel_pos;
        std::copy(vpos.begin(), vpos.end(), pos.begin());
    }

    void update_fields(const VecSlice& vel, const RealSlice& density) override
    {
        std::copy(vel.begin(), vel.end(), m_data.grid().vel.begin());
        std::copy(
            density.begin(), density.end(), m_data.grid().density.begin());
        ops::UpdateVelOp<ActTrait, SrcTrait>()(m_data);
    }

    void compute_forces() override
    {
        ops::ComputeForceOp<ActTrait, SrcTrait>()(m_data);
        m_src_op.setup_op();
    }

    void compute_source_term(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Geometry& geom) override
    {
        m_src_op(lev, mfi, geom);
    }

    void prepare_outputs(const std::string& out_dir) override
    {
        m_out_op.prepare_outputs(out_dir);
    }

    void write_outputs() override { m_out_op.write_outputs(); }

    void init_actuator_source() override
    {
        ops::InitDataOp<ActTrait, SrcTrait>()(m_data);
        m_src_op.initialize();
    }
};

template <typename ActTrait, typename SrcTrait>
void ActModel<ActTrait, SrcTrait>::determine_root_proc(
    amrex::Vector<int>& act_proc_count)
{
    ops::determine_root_proc<ActTrait>(m_data, act_proc_count);
    {
        // Sanity checks
        const auto& linfo = m_data.info();
        const auto& plist = linfo.procs;
        AMREX_ALWAYS_ASSERT(linfo.root_proc > -1);
        AMREX_ALWAYS_ASSERT(plist.find(linfo.root_proc) != plist.end());
    }
}

template <typename ActTrait, typename SrcTrait>
int ActModel<ActTrait, SrcTrait>::num_velocity_points() const
{
    const auto& linfo = m_data.info();
    auto& grid = m_data.grid();

    return (linfo.sample_vel_in_proc) ? grid.vel.size() : 0;
}

} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATORMODEL_H */
