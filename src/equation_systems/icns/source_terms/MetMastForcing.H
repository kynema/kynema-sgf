#ifndef METMASTFORCING_H
#define METMASTFORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"

namespace kynema_sgf::pde::icns {

/** Relax the velocity towards met-mast and lidar measurements
 *  \ingroup icns_src
 *
 *  Each measurement station is a horizontal location and a profile of one or
 *  more heights above the terrain. A met-mast point is a one-level profile and
 *  a lidar is a multi-level profile with the standard deviation of the wind.
 *  Station \f$i\f$ has the weight
 *
 *  \f[ w_i = \exp\left(-\frac{1}{4}\left[\frac{r_i^2}{R_h^2} +
 *  \frac{d_i^2}{R_z^2}\right]\right) \f]
 *
 *  with \f$r_i\f$ the horizontal distance to the station and \f$d_i\f$ the
 *  vertical distance outside the measured height range (zero inside it). The
 *  target is the measured velocity \f$\bar{u}_i\f$ widened to a band of
 *  \f$\pm\alpha \sigma_i\f$, so that only the part of the velocity outside
 *  the measured variability is forced:
 *
 *  \f[ t_i = \min(\max(u, \bar{u}_i - \alpha\sigma_i), \bar{u}_i +
 *  \alpha\sigma_i) \f]
 *
 *  The stations are combined into one relaxation
 *
 *  \f[ S = -\min\left(\frac{W}{\tau}, \frac{1}{\Delta t}\right) (u - \bar{t}),
 *  \quad W = \min\left(1, \sum_i w_i\right), \quad \bar{t} = \frac{\sum_i w_i
 *  t_i}{\sum_i w_i} \f]
 *
 *  The heights are above the local terrain when TerrainDrag is active, and no
 *  forcing is applied in the cells inside the terrain.
 */
class MetMastForcing : public MomentumSource::Register<MetMastForcing>
{
public:
    static std::string identifier() { return "MetMastForcing"; }

    explicit MetMastForcing(const CFDSim& sim);

    ~MetMastForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

    //! Number of measurement stations (met-mast points and lidars)
    [[nodiscard]] int num_stations() const
    {
        return static_cast<int>(m_station_x.size());
    }

private:
    //! Read met-mast points, one per line: x y z u v w T
    void read_point_file(const std::string& fname);

    //! Read a lidar profile: x y on the first line, then z u v w su sv sw
    void read_profile_file(const std::string& fname);

    void add_level(
        amrex::Real z,
        const amrex::Array<amrex::Real, AMREX_SPACEDIM>& vel,
        const amrex::Array<amrex::Real, AMREX_SPACEDIM>& sigma);

    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const CFDSim& m_sim;

    //! Station locations
    amrex::Vector<amrex::Real> m_station_x;
    amrex::Vector<amrex::Real> m_station_y;
    //! Levels of station i are [m_level_offset[i], m_level_offset[i+1])
    amrex::Vector<int> m_level_offset{0};
    //! Height above terrain of each level
    amrex::Vector<amrex::Real> m_level_z;
    //! Velocity and its standard deviation, AMREX_SPACEDIM values per level
    amrex::Vector<amrex::Real> m_level_vel;
    amrex::Vector<amrex::Real> m_level_sigma;

    amrex::Gpu::DeviceVector<amrex::Real> m_station_x_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_station_y_d;
    amrex::Gpu::DeviceVector<int> m_level_offset_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_level_z_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_level_vel_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_level_sigma_d;

    //! Relaxation time scale
    amrex::Real m_timescale{30};
    //! Horizontal radius of the forcing footprint
    amrex::Real m_horizontal_radius{500};
    //! Vertical radius of the forcing footprint
    amrex::Real m_vertical_radius{25};
    //! Cutoff on the normalized squared distance
    amrex::Real m_damping_radius{1400};
    //! Half-width of the tolerance band in standard deviations
    amrex::Real m_sigma_factor{1};
};

} // namespace kynema_sgf::pde::icns
#endif
