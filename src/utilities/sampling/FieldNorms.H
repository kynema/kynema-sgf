#ifndef FIELDNORMS_H
#define FIELDNORMS_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"

namespace amr_wind::field_norms {

/** field norm object
 *  \ingroup sampling
 *
 *  A concrete implementation of the post-processing interface that deals with
 *  norms of fields. FIXME: need to update this
 *
 */
class FieldNorms : public PostProcessBase::Register<FieldNorms>
{
public:
    static std::string identifier() { return "FieldNorms"; }

    FieldNorms(CFDSim& /*sim*/, std::string /*label*/);

    ~FieldNorms() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs and create the different data probe instances
    void initialize() override;

    void post_advance_work() override {};

    void output_actions() override;

    void post_regrid_actions() override {}

    static amrex::Real get_norm(
        const amr_wind::Field& field,
        int comp,
        int ncomp,
        int norm_type,
        bool use_mask);

    [[nodiscard]] const amrex::Vector<std::string>& var_names() const
    {
        return m_var_names;
    }

    [[nodiscard]] const amrex::Vector<amrex::Real>& field_norms() const
    {
        return m_fnorms;
    }

private:
    //! prepare ASCII file and directory
    virtual void prepare_ascii_file();

    //! Output sampled data in ASCII format
    virtual void write_ascii();

    //!  process all fields and calculate norms
    void process_field_norms();

    //! Reference to the CFD sim
    CFDSim& m_sim;

    //! List of variable names for output
    amrex::Vector<std::string> m_var_names;

    //! List holding norms for all fields and their components
    amrex::Vector<amrex::Real> m_fnorms;

    /** Name of this sampling object.
     *
     *  The label is used to read user inputs from file and is also used for
     *  naming files directories depending on the output format.
     */
    const std::string m_label;

    //! filename for ASCII output
    std::string m_out_fname;

    //! width in ASCII output
    int m_width{18};

    //! precision in ASCII output
    int m_precision{10};

    //! Type of norm
    int m_norm_type{2};

    //! Mask redundant grids (finest available grid only)
    bool m_use_mask{true};

    //! Consider vector magnitude, not separate vector components
    bool m_use_vector_magnitude{false};
};

} // namespace amr_wind::field_norms

#endif /* FIELDNORMS_H */
