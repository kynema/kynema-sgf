#ifndef KINETICENERGY_H
#define KINETICENERGY_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::kinetic_energy {

/** kinetic energy object
 *  \ingroup sampling
 *
 *  A concrete implementation of the post-processing interface that deals with
 *  integrating total kinetic energy
 *
 */
class KineticEnergy : public PostProcessBase::Register<KineticEnergy>
{
public:
    static std::string identifier() { return "KineticEnergy"; }

    KineticEnergy(CFDSim& /*sim*/, std::string /*label*/);

    ~KineticEnergy() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and create the different data probe instances
    void initialize() override;

    void post_advance_work() override {};

    void output_actions() override;

    void post_regrid_actions() override {}

    //! calculate the L2 norm of a given field and component
    amrex::Real calculate_kinetic_energy();

private:
    //! prepare ASCII file and directory
    void prepare_ascii_file();

    //! Output sampled data in ASCII format
    void write_ascii();

    //! store the total kinetic energy
    amrex::Real m_total_kinetic_energy{0.0_rt};

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
    //! reference to density
    const Field& m_density;

    //! filename for ASCII output
    std::string m_out_fname;

    //! width in ASCII output
    int m_width{22};

    //! precision in ASCII output
    int m_precision{12};
};

} // namespace amr_wind::kinetic_energy

#endif /* KINETICENERGY_H */
