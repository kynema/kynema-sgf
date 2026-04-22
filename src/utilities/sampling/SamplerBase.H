#ifndef SAMPLERBASE_H
#define SAMPLERBASE_H

#include <AMReX_RealVect.H>
#include "amr-wind/core/Factory.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class CFDSim;

namespace sampling {

struct SampleLocType
{
    using LocType = amrex::Vector<amrex::RealVect>;
    using IdType = amrex::Vector<amrex::Long>;

    void push_back(const amrex::RealVect& loc, const amrex::Long id)
    {
        m_locations.push_back(loc);
        m_ids.push_back(id);
    }

    [[nodiscard]] const LocType& locations() const
    {
        AMREX_ALWAYS_ASSERT(m_locations.size() == m_ids.size());
        return m_locations;
    }

    LocType& locations()
    {
        AMREX_ALWAYS_ASSERT(m_locations.size() == m_ids.size());
        return m_locations;
    }

    [[nodiscard]] const IdType& ids() const
    {
        AMREX_ALWAYS_ASSERT(m_locations.size() == m_ids.size());
        return m_ids;
    }

    LocType m_locations;
    IdType m_ids;
};

/** Abstract representation of data probes to sample flow data
 *  \ingroup sampling
 *
 *  This class defines the interface by which the positions of the data-probes
 *  are defined and access from the Sampling class. All data probe definitions
 *  must be sub-classes of this base class.
 */
class SamplerBase : public Factory<SamplerBase, CFDSim&>
{
public:
    static std::string base_identifier() { return "SamplerBase"; }

    static constexpr amrex::Real bounds_tol = constants::TIGHT_TOL;

    ~SamplerBase() override = default;

    //! Name used to refer to this sampler (e.g., myline1)
    [[nodiscard]] virtual std::string label() const = 0;
    virtual std::string& label() = 0;

    //! Class name of this sampler (e.g., LineSampler)
    [[nodiscard]] virtual std::string sampletype() const = 0;

    //! Unique numeric identifier used to track particles belonging to this
    //! instance
    [[nodiscard]] virtual int id() const = 0;
    virtual int& id() = 0;

    //! Total number of probes that belong to this sampling instance
    [[nodiscard]] virtual long num_points() const = 0;

    //! Total number of probes to be output
    //! Note: Same as num_points() for sampling without data reduction
    [[nodiscard]] virtual long num_output_points() const = 0;

    //! Read inputs and perform initialization actions
    virtual void initialize(const std::string& key) = 0;

    //! Populate the vector with coordinates of the sampling locations
    virtual void sampling_locations(SampleLocType&) const = 0;

    //! Populate the vector with coordinates of the sampling locations inside
    //! a box
    virtual void
    sampling_locations(SampleLocType&, const amrex::Box&) const = 0;

    //! Check and fix the bounds of the sampler so the probes are in the domain
    virtual void check_bounds() = 0;

    //! Populate the vector with coordinates of the output locations
    virtual void output_locations(SampleLocType&) const = 0;

    //! Update the sampling locations
    virtual bool update_sampling_locations() { return false; }

    //! Run actions after sample (useful in interpolated subsampling)
    virtual void post_sample_actions() {}

    //! Run actions after regrid (important when using field-based quantities)
    virtual void post_regrid_actions() {}

    //! Run specific output for the sampler
    virtual bool output_netcdf_field(
        const std::vector<amrex::Real>& /*unused*/,
        ncutils::NCGroup& /*unused*/,
        const size_t /*unused*/)
    {
        return true;
    }

    //! Run data modification for specific sampler?
    virtual bool do_data_modification() { return false; }
    virtual bool do_convert_velocity_los() { return false; }
    virtual bool do_subsampling_interp() { return false; }

    //! Sample buffer modification instructions
    virtual std::vector<amrex::Real> modify_sample_data(
        const std::vector<amrex::Real>& sampledata,
        const std::string& /*unused*/)
    {
        return sampledata;
    }

    virtual void calc_lineofsight_velocity(
        const std::vector<std::vector<amrex::Real>>& /*unused*/,
        const int /*unused*/)
    {}

    //! Populate metadata in the NetCDF file
    virtual void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const
    {}
    virtual void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const
    {}
    virtual void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
    {}
};

} // namespace sampling
} // namespace amr_wind

#endif /* SAMPLERBASE_H */
