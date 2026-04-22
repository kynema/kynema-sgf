#ifndef BODYFORCE_H
#define BODYFORCE_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "AMReX_GpuContainers.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include <cstddef>

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Apply a user-defined body force as a source term
 *  \ingroup icns_src
 *
 *  The implementation can apply either a constant source term, or a sinusoidal
 *  term as a function of time.
 */
class BodyForce : public MomentumSource::Register<BodyForce>
{
public:
    static std::string identifier() { return "BodyForce"; }

    explicit BodyForce(const CFDSim& /*sim*/);

    ~BodyForce() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void read_bforce_profile(const std::string& filename);
    void read_bforce_timetable(const std::string& filename);

private:
    //! Time
    const SimTime& m_time;

    const amrex::AmrCore& m_mesh;

    //! Forcing source term (pressure gradient)
    amrex::Vector<amrex::Real> m_body_force{0.0_rt, 0.0_rt, 0.0_rt};

    //! Body Force Type
    std::string m_type{"uniform_constant"};
    //! Uniform time table file
    std::string m_utt_file;

    //! Angular frequency used in the oscillatory forcing
    amrex::Real m_omega{0.0_rt};
    std::string m_bforce_file;
    size_t m_bforce_profile_nhts;

    amrex::Gpu::DeviceVector<amrex::Real> m_prof_x;
    amrex::Gpu::DeviceVector<amrex::Real> m_prof_y;
    amrex::Gpu::DeviceVector<amrex::Real> m_ht;

    //! Vectors for storing uniform_timetable inputs
    amrex::Vector<amrex::Real> m_time_table;
    amrex::Vector<amrex::Real> m_fx_table;
    amrex::Vector<amrex::Real> m_fy_table;
    amrex::Vector<amrex::Real> m_fz_table;
};

} // namespace amr_wind::pde::icns

#endif /* BODYFORCE_H */
