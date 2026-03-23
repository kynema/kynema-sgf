#ifndef ABLMESOSCALEFORCING_H
#define ABLMESOSCALEFORCING_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "AMReX_Vector.H"
#include "amr-wind/CFDSim.H"

using namespace amrex::literals;

namespace amr_wind {

class ABLMesoscaleForcing
{

public:
    ABLMesoscaleForcing(const CFDSim& sim, const std::string& identifier);

    void indirect_forcing_init();

    static void invert_mat(
        const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
        amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im);

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

    int m_idx_time{0};
    int m_axis{2};

    // inputs
    std::string m_forcing_scheme;
    amrex::Real m_gain_coeff{
        0.2_rt};            // relaxation factor for all forcing schemes
    bool m_debug{false};    // print out profiles for debugging
    bool m_tendency{false}; // check if tendency forcing

    // planar averaging heights
    amrex::Vector<amrex::Real> m_zht;
    int m_nht{0}; // set by derived ABLMesoForcing[Mom|Temp] class

    //
    // regression vars for indirect forcing
    //
    // int m_ind_polyOrder;
    int m_norm_zmax = 0;
    amrex::Real m_scaleFact =
        1.0e-3_rt; // to improve conditioning of regression matrix

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> m_im_zTz;

    // - regression weighting profile inputs (these may be constant or the
    // instantaneous value, if forcing_transition is specified)
    std::string m_forcing_transition;
    amrex::Real m_transition_height;
    bool m_update_transition_height{false};
    amrex::Real m_transition_thickness;

    amrex::Vector<amrex::Real> m_weighting_heights;
    amrex::Vector<amrex::Real> m_weighting_values;
    bool m_user_specified_weighting{false};

    // - W(z) at planar averaging heights, for regression
    amrex::Vector<amrex::Real> m_W;

    // - blending heights (constant over time if the transition height is also
    // constant)
    amrex::Vector<amrex::Real> m_blending_heights;
    amrex::Vector<amrex::Real> m_blending_values;

    // - weights at planar averaging heights, for blending indirect-direct
    // forcing profiles
    amrex::Vector<amrex::Real> m_blend;

    //
    // weighting profile manipulation
    //
    void set_transition_weighting(); // will override weighting profile
    void
    update_weights(); // interpolate current profile to planar averaging heights

    //
    // forcing profile manipulation
    //
    bool have_forcing_transition()
    {
        return (amrex::toLower(m_forcing_transition) != "none");
    }

    bool forcing_to_constant()
    {
        // - check if m_forcing_transition ends with "Constant"
        if (m_forcing_transition.length() < 8) {
            return false;
        }
        return (m_forcing_transition.ends_with("Constant"));
    }

    // - slope varies linearly from the actual gradient of the forcing profile
    // at m_transition_height
    //   to 0 at (m_transition_height + m_transition_thicnkess)
    void constant_forcing_transition(amrex::Vector<amrex::Real>& error);

    // - blend between a lower and upper forcing profile; the blending fraction
    // is assumed to be
    //   equal to W(z), but this is not a requirement; W=1 and W=0 correspond to
    //   the lower and upper profiles, respectively
    void blend_forcings(
        const amrex::Vector<amrex::Real>& lower,
        const amrex::Vector<amrex::Real>& upper,
        amrex::Vector<amrex::Real>& error);
};

} // namespace amr_wind

#endif
