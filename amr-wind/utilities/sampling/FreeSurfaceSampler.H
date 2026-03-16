#ifndef FREESURFACESAMPLER_H
#define FREESURFACESAMPLER_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind::sampling {

class FreeSurfaceSampler : public SamplerBase::Register<FreeSurfaceSampler>
{
public:
    static std::string identifier() { return "FreeSurfaceSampler"; }

    FreeSurfaceSampler(CFDSim& /*sim*/);

    ~FreeSurfaceSampler() override;

    //! Read user inputs and create mesh-based information for locating
    //! interface
    void initialize(const std::string& key) override;

    //! Check and fix the bounds of the sampler so the probes are in the domain
    void check_bounds() override;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType& /*sample_locs*/) const override;

    //! Populate the vector with coordinates of the sampling locations inside
    //! a box
    void sampling_locations(
        SampleLocType& /*sample_locs*/,
        const amrex::Box& /*box*/) const override;

    void output_locations(SampleLocType& sample_locs) const override
    {
        sampling_locations(sample_locs);
    }

    //! Find heights associated with 2D sample locations
    bool update_sampling_locations() override;

    //! Redo some of the initialization work when the grid changes
    void post_regrid_actions() override;

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/, size_t /*unused*/) const override;

    //! Name of this sampling object
    [[nodiscard]] std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this sampling object
    [[nodiscard]] std::string sampletype() const override
    {
        return identifier();
    }

    //! Unique identifier for this set of probe locations
    [[nodiscard]] int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Output functions for private variables
    //! Number of 2D grid points
    [[nodiscard]] int num_gridpoints() const { return m_npts; }
    //! Number of instances
    [[nodiscard]] int num_instances() const { return m_ninst; }
    //! Number of points, total
    [[nodiscard]] long num_points() const override
    {
        return static_cast<long>(m_npts) * m_ninst;
    }
    [[nodiscard]] long num_output_points() const override
    {
        return static_cast<long>(m_npts) * m_ninst;
    }
    //! Locations
    [[nodiscard]] amrex::Vector<amrex::Array<amrex::Real, 2>>
    grid_locations() const
    {
        return m_grid_locs;
    }
    //! Outputs (heights)
    [[nodiscard]] amrex::Vector<amrex::Real> heights() const { return m_out; }

private:
    CFDSim& m_sim;

    //! reference to VOF
    const Field& m_vof;

    //! Number or points on 2D grid in each direction
    amrex::Vector<int> m_npts_dir;
    int m_npts{0};
    //! Number of instances (possible sampling points per location, like in the
    //! case of a breaking wave)
    int m_ninst{1};

    //! Coordinate direction to search along, default is z
    //! (this also determines the meaning of start and end points)
    int m_coorddir{2};
    //! Grid coordinates, determined as a function of m_coorddir
    int m_gc0 = 0;
    int m_gc1 = 1;

    //! Parameters to set up plane
    amrex::Vector<amrex::Real> m_start, m_end;
    //! Locations of points in 2D grid
    amrex::Vector<amrex::Array<amrex::Real, 2>> m_grid_locs;
    //! Output coordinate
    amrex::Vector<amrex::Real> m_out;

    //! Flag to use linear interpolation for surface finding in part of domain
    bool m_use_linear{false};
    //! Parameter for optional linear interpolation
    amrex::Real m_lx_linear{0.};

    //! Max number of sample points found in a single cell
    int m_ncomp{1};
    //! Max number of sample points allowed in a single cell
    int m_ncmax{8};

    std::string m_label;
    int m_id{-1};
};

} // namespace amr_wind::sampling

#endif /* FREESURFACESAMPLER_H */
