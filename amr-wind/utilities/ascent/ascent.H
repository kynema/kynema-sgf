#ifndef ASCENT_INT_H
#define ASCENT_INT_H

#include "amr-wind/utilities/PostProcessing.H"

/**
 * Ascent In-situ Integration
 */

namespace amr_wind {

class Field;

namespace ascent_int {

class AscentPostProcess : public PostProcessBase::Register<AscentPostProcess>
{
public:
    static const std::string identifier() { return "Ascent"; }

    AscentPostProcess(CFDSim&, const std::string&);

    virtual ~AscentPostProcess();

    //! Perform actions before mesh is created
    void pre_init_actions() override;

    //! Read user inputs and create the different data probe instances
    void initialize() override;

    void post_advance_work() override {};

    void output_actions() override;

    //! Actions to perform post regrid e.g. redistribute particles
    void post_regrid_actions() override;

protected:
private:
    CFDSim& m_sim;
    std::string m_label;

    amrex::Vector<std::string> m_var_names;
    amrex::Vector<Field*> m_fields;
};

} // namespace ascent_int
} // namespace amr_wind

#endif
