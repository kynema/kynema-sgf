#ifndef MODULATEDPOWERLAW_H
#define MODULATEDPOWERLAW_H

#include "src/core/Field.H"
#include "src/CFDSim.H"
#include "src/boundary_conditions/field_boundary_fill/FieldBoundary.H"
#include "AMReX_Gpu.H"
#include "src/utilities/trig_ops.H"
#include "src/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

/** Interface for ABL modulated power law
 *  \ingroup we_abl
 *
 *  This class sets velocity and temperature at mass inflow bc's
 *
 */
class ModulatedPowerLaw : public FieldBoundary::Register<ModulatedPowerLaw>
{

public:
    static std::string identifier() { return "ModulatedPowerLaw"; }

    explicit ModulatedPowerLaw(CFDSim& /*sim*/);

    void post_init_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

    void set_velocity(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab,
        int dcomp = 0,
        int orig_comp = 0) const override;

    void set_temperature(
        int lev,
        amrex::Real time,
        const Field& fld,
        amrex::MultiFab& mfab) const;

private:
    const CFDSim& m_sim;
    const kynema_sgf::SimTime& m_time;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    Field& m_temperature;

    amrex::Real m_zoffset{0.0_rt};
    amrex::Real m_zref{90};
    amrex::Real m_shear_exp{0.1_rt};
    amrex::Real m_umax_factor{1.2_rt};
    amrex::Real m_bulk_velocity{15.0_rt};
    amrex::Real m_shearlayer_height{600.0_rt};
    amrex::Real m_shearlayer_smear_thickness{30.0_rt};
    amrex::Real m_wind_speed{8.0_rt};
    amrex::Real m_wind_direction{270.0_rt};
    amrex::Vector<amrex::Real> m_uvec{8.0_rt, 0.0_rt, 0.0_rt};

    amrex::Real m_start_time{0.0_rt};
    amrex::Real m_stop_time{std::numeric_limits<amrex::Real>::max()};
    amrex::Real m_degrees_per_sec{0.02_rt};

    ///@{
    //! Lookup table containing temperature variation w.r.t. height above
    //! terrain
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;
    ///@}

    // Device copies of the above arrays
    amrex::Gpu::DeviceVector<amrex::Real> m_thht_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_thvv_d;

    amrex::Real m_delta_t{0.8_rt};
    amrex::Real m_theta_cutoff_height{250.0_rt};
    amrex::Real m_theta_gauss_mean{0.0_rt};
    amrex::Real m_theta_gauss_var{1.0_rt};
};

} // namespace kynema_sgf

#endif /* MODULATEDPOWERLAW_H */
