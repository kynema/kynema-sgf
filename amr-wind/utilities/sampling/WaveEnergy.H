#ifndef WAVEENERGY_H
#define WAVEENERGY_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::wave_energy {

/** wave energy object
 *  \ingroup sampling
 *
 *  A concrete implementation of the post-processing interface that deals with
 *  integrating total wave energy. This routine calculates the kinetic energy
 *  (KE) and potential energy (PE) according to the following definitions. It
 *  calculates the correct potential energy offset, assuming an equilibrium
 *  water level of 0 unless specified otherwise by an input parameter. It also
 *  assumes the traditional coordinate system, where interface height/depth is
 *  in z.
 *
 *  Integral definitions (not output)
 *  KE_i = 0.5_rt * (Integral over liquid phase) \\rho * u^2 dV
 *  PE_i = (Integral over liquid phase) \\rho * g * z_liq dV
 *  (offset for PE would be 1/2 * Lx * Ly * \\rho_liq * g * depth^2)
 *
 *  Normalized definitions (which are output) in form of avg kinematic energy
 *  KE = KE_i / (d * \\rho_liq * Lx * Ly)
 *  PE = PE_i / (d * \\rho_liq * Lx * Ly) + 1/2 * g * depth
 */
class WaveEnergy : public PostProcessBase::Register<WaveEnergy>
{
public:
    static std::string identifier() { return "WaveEnergy"; }

    WaveEnergy(CFDSim& /*sim*/, std::string /*label*/);

    ~WaveEnergy() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and get information needed for calculations
    void initialize() override;

    void post_advance_work() override {};

    //! Integrate energy components and output to file
    void output_actions() override;

    void post_regrid_actions() override {}

    //! Calculate the sum of stated energy in liquid phase
    amrex::Real calculate_kinetic_energy();
    amrex::Real calculate_potential_energy();

    //! Output private variables that store energy measurements
    void wave_energy(amrex::Real& ke, amrex::Real& pe) const
    {
        ke = m_wave_kinetic_energy;
        pe = m_wave_potential_energy;
    };

private:
    //! prepare ASCII file and directory
    virtual void prepare_ascii_file();

    //! Output sampled data in ASCII format
    virtual void write_ascii();

    //! store the total wave energy
    amrex::Real m_wave_kinetic_energy{0.0_rt};
    amrex::Real m_wave_potential_energy{0.0_rt};

    //! Reference to the CFD sim
    CFDSim& m_sim;

    /** Name of this sampling object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories depending on the output format.
     */
    const std::string m_label;

    //! reference to velocity
    const Field& m_velocity;
    //! reference to vof
    const Field& m_vof;

    //! gravity vector
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    //! offset for potential energy calculation
    amrex::Real m_pe_off = 0.0_rt;
    //! volumetric scaling for each energy
    amrex::Real m_escl = 1.0_rt;

    //! filename for ASCII output
    std::string m_out_fname;

    //! width in ASCII output
    int m_width{22};

    //! precision in ASCII output
    int m_precision{12};
};

} // namespace amr_wind::wave_energy

#endif /* WAVEENERGY_H */
