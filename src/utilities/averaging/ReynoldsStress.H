#ifndef REYNOLDSSTRESS_H
#define REYNOLDSSTRESS_H

#include "amr-wind/utilities/averaging/TimeAveraging.H"

namespace amr_wind::averaging {

/** Compute correlation from 2 CFD fields
 *  The formula used to compute the correlations between two fields:
 *  A = < A > + a and B = < B > + b is < ab > = < AB > -< A >< B >
 *  where A and B are the mean values and a and b are the fluctuations
 *
 */
class ReynoldsStress : public FieldTimeAverage::Register<ReynoldsStress>
{
public:
    static std::string identifier() { return "ReynoldsStress"; }

    ReynoldsStress(
        CFDSim& /*sim*/, const std::string& avgname, const std::string& fname);

    /** Update field averaging at a given timestep
     *
     *  \param time SimTime instance
     *  \param filter_width Time-averaging window specified by user
     *  \param avg_time_interval Time interval at which averaging is performed
     *  \param elapsed_time Time elapsed since averaging was initiated
     */
    void operator()(
        const SimTime& /*time*/,
        amrex::Real /*filter_width*/,
        amrex::Real /*avg_time_interval*/,
        amrex::Real /*elapsed_time*/) override;

    const std::string& average_field_name() override;

private:
    //! Generate the Reynolds stress field name based on the field name and
    //! label
    static std::string
    stress_name(const std::string& fname, const std::string& avgname)
    {
        return fname + "_" + "stress" + "_" + avgname;
    }

    static std::string
    re_stress_name(const std::string& fname, const std::string& avgname)
    {
        return fname + "_" + "reynolds_stress" + "_" + avgname;
    }

    //! Fluctuating field
    const Field& m_field;

    //! Reynolds averaged field
    const Field& m_average;

    //! The stresses < AB >
    Field& m_stress;

    //! The reynolds stresses < ab > = < AB > - < A > < B >
    Field& m_re_stress;
};

} // namespace amr_wind::averaging

#endif /* REYNOLDSSTRESS_H */
