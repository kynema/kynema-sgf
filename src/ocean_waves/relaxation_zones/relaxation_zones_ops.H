#ifndef RELAX_ZONES_OPS_H
#define RELAX_ZONES_OPS_H

#include "src/ocean_waves/relaxation_zones/RelaxationZones.H"
#include "src/ocean_waves/OceanWavesTypes.H"
#include "src/ocean_waves/OceanWavesOps.H"
#include "src/ocean_waves/utils/wave_utils_K.H"
#include "src/core/MultiParser.H"
#include "src/fvm/gradient.H"

#include "AMReX_MultiFabUtil.H"

namespace kynema_sgf::ocean_waves {
namespace relaxation_zones {

/** Parse common data for modeling relaxation zones
 */
void read_inputs(
    RelaxZonesBaseData& /*wdata*/,
    OceanWavesInfo& /*unused*/,
    const ::kynema_sgf::utils::MultiParser& /*pp*/);

/** Initialize core data structures when modeling the IB geometry
 */
void init_data_structures(CFDSim&);

/** Harmonize the target wave solution with numerical beach */
void modify_target_fields_for_beach(
    CFDSim& sim, const RelaxZonesBaseData& wdata);

/** Set the free_surface height inside the relaxation zones
 */
void apply_relaxation_zones(CFDSim& sim, const RelaxZonesBaseData& wdata);

void update_target_vof(CFDSim& sim);

void prepare_netcdf_file(
    const std::string& /*ncfile*/,
    const RelaxZonesBaseData& /*meta*/,
    const OceanWavesInfo& /*info*/);

void write_netcdf(
    const std::string& /*ncfile*/,
    const RelaxZonesBaseData& /*meta*/,
    const OceanWavesInfo& /*info*/,
    amrex::Real /*time*/);

} // namespace relaxation_zones

namespace ops {

template <typename WaveTheoryTrait>
struct UseDefaultOp
{
    static constexpr bool update_relaxation_zones = true;
    static constexpr bool apply_relaxation_zones = true;
};

template <typename WaveTheoryTrait>
struct UpdateTargetFieldsOp<
    WaveTheoryTrait,
    typename std::enable_if_t<
        std::is_base_of_v<RelaxZonesType, WaveTheoryTrait>>>
{
    void operator()(
        typename WaveTheoryTrait::DataType& /*unused*/,
        const amrex::Real /*unused*/)
    {}
};

template <typename WaveTheoryTrait>
struct ApplyRelaxZonesOp<
    WaveTheoryTrait,
    typename std::enable_if_t<
        std::is_base_of_v<RelaxZonesType, WaveTheoryTrait>>>
{
    void operator()(typename WaveTheoryTrait::DataType& data)
    {
        BL_PROFILE("amr-wind::ocean_waves::OceanWaves::update_velocities");

        const auto& wdata = data.meta();
        auto& sim = data.sim();

        relaxation_zones::apply_relaxation_zones(sim, wdata);
    }
};

template <typename WaveTheoryTrait>
struct UpdateTargetFieldsSecondaryOp<
    WaveTheoryTrait,
    typename std::enable_if_t<
        std::is_base_of_v<RelaxZonesType, WaveTheoryTrait>>>
{
    void operator()(typename WaveTheoryTrait::DataType& data)
    {
        BL_PROFILE("amr-wind::ocean_waves::OceanWaves::update_target_vof");

        const auto& wdata = data.meta();
        auto& sim = data.sim();

        if (wdata.has_beach && sim.repo().field_exists("vof")) {
            relaxation_zones::modify_target_fields_for_beach(sim, wdata);
        }

        relaxation_zones::update_target_vof(sim);
    }
};

template <typename WaveTheoryTrait>
struct ProcessOutputsOp<
    WaveTheoryTrait,
    typename std::enable_if_t<
        std::is_base_of_v<RelaxZonesType, WaveTheoryTrait>>>
{
private:
    typename WaveTheoryTrait::DataType& m_data;

    std::string m_out_dir;
    std::string m_nc_filename;
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename WaveTheoryTrait::DataType& data)
        : m_data(data)
    {}

    void read_io_options(const ::kynema_sgf::utils::MultiParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        relaxation_zones::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        relaxation_zones::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), time.new_time());
    }
};

} // namespace ops

} // namespace kynema_sgf::ocean_waves

#endif /* RELAX_ZONES_OPS_H */
