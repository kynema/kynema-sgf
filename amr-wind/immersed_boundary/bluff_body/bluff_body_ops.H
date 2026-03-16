#ifndef BLUFF_BODY_OPS_H
#define BLUFF_BODY_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"
#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/fvm/gradient.H"

#include "AMReX_MultiFabUtil.H"

namespace amr_wind::ib {
namespace bluff_body {

/** Parse common data for modeling non-blade bluff bodies via an immersed
 * boundary
 */
void read_inputs(
    BluffBodyBaseData& /*wdata*/,
    IBInfo& /*unused*/,
    const ::amr_wind::utils::MultiParser& /*pp*/);

/** Initialize core data structures when modeling the IB geometry
 */
void init_data_structures(CFDSim&);

/** Set the velocity inside the IB based on a manufactured solution
 */
void apply_mms_vel(CFDSim& /*sim*/);

/** Set the velocity inside the IB based on a dirichlet BC
 */
void apply_dirichlet_vel(
    CFDSim& /*sim*/, const amrex::Vector<amrex::Real>& vel_bc);

void prepare_netcdf_file(
    const std::string& /*ncfile*/,
    const BluffBodyBaseData& /*meta*/,
    const IBInfo& /*info*/);

void write_netcdf(
    const std::string& /*ncfile*/,
    const BluffBodyBaseData& /*meta*/,
    const IBInfo& /*info*/,
    amrex::Real /*time*/);

} // namespace bluff_body

namespace ops {

template <typename GeomTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename GeomTrait>
struct UpdatePosOp<
    GeomTrait,
    typename std::enable_if_t<std::is_base_of_v<BluffBodyType, GeomTrait>>>
{
    void operator()(typename GeomTrait::DataType& /*unused*/) {}
};

template <typename GeomTrait>
struct UpdateVelOp<
    GeomTrait,
    typename std::enable_if_t<std::is_base_of_v<BluffBodyType, GeomTrait>>>
{
    void operator()(typename GeomTrait::DataType& data)
    {
        BL_PROFILE("amr-wind::ib::IB::update_velocities");

        const auto& wdata = data.meta();
        auto& sim = data.sim();

        if (wdata.is_mms) {
            bluff_body::apply_mms_vel(sim);
        } else {
            bluff_body::apply_dirichlet_vel(sim, wdata.vel_bc);
        }
    }
};

template <typename GeomTrait>
struct ComputeForceOp<
    GeomTrait,
    typename std::enable_if_t<std::is_base_of_v<BluffBodyType, GeomTrait>>>
{
    void operator()(typename GeomTrait::DataType& /*unused*/) {}
};

template <typename GeomTrait>
struct ProcessOutputsOp<
    GeomTrait,
    typename std::enable_if_t<std::is_base_of_v<BluffBodyType, GeomTrait>>>
{
private:
    typename GeomTrait::DataType& m_data;

    std::string m_out_dir;
    std::string m_nc_filename;
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename GeomTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const ::amr_wind::utils::MultiParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        bluff_body::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        bluff_body::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), time.new_time());
    }
};

} // namespace ops

} // namespace amr_wind::ib

#endif /* BLUFF_BODY_OPS_H */
