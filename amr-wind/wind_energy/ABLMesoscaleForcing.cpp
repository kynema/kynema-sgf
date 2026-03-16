#include "amr-wind/wind_energy/ABLMesoscaleForcing.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

// WORKAROUND
#include <fstream>

using namespace amrex::literals;

namespace amr_wind {

ABLMesoscaleForcing::ABLMesoscaleForcing(
    const CFDSim& sim, const std::string& identifier)
    : m_time(sim.time()), m_mesh(sim.mesh())
{
    amrex::Print() << "Constructing " << identifier << " object" << '\n';

    amrex::ParmParse pp(identifier);
    pp.query("forcing_scheme", m_forcing_scheme);
    pp.query("control_gain", m_gain_coeff);
    pp.query("debug", m_debug);
    amrex::Print() << "  forcing_scheme : " << m_forcing_scheme << '\n';
    amrex::Print() << "  control_gain   : " << m_gain_coeff << '\n';

    if (pp.query("forcing_transition", m_forcing_transition) == 0) {
        amrex::Print() << "  using full profile assimilation by default"
                       << '\n';
        m_forcing_transition = "none";
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {

        // check for user-specified weighting profile
        if (pp.queryarr("weighting_heights", m_weighting_heights) == 1) {
            pp.getarr("weighting_values", m_weighting_values);
            amrex::Print() << "  given weighting profile:" << '\n';
            for (int i = 0; i < m_weighting_heights.size(); ++i) {
                amrex::Print() << "  " << m_weighting_heights[i] << " "
                               << m_weighting_values[i] << '\n';
            }
            AMREX_ALWAYS_ASSERT(
                m_weighting_heights.size() == m_weighting_values.size());
            m_user_specified_weighting = true;
        }

        if (!have_forcing_transition()) {
            if (!m_user_specified_weighting) {
                // default is to have uniform weighting throughout
                amrex::Print()
                    << "  using default regression weighting" << '\n';
                amrex::Real zmin = m_mesh.Geom(0).ProbLo(m_axis);
                amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
                m_weighting_heights = {zmin, zmax};
                m_weighting_values = {1.0_rt, 1.0_rt};
            }
        } else {
            // weightings will be automatically set based on forcing transition
            pp.get(
                "transition_thickness",
                m_transition_thickness); // constant, required
            if (pp.query("constant_transition_height", m_transition_height) ==
                1) {
                amrex::Print()
                    << "  fixed transition layer between "
                    << m_transition_height << " and "
                    << m_transition_height + m_transition_thickness << '\n';
                // set weighting profile
                set_transition_weighting();
            } else {
                // expect to read transition_height history in netCDF input file
                // and the weighting profile will need to be updated every step
                m_update_transition_height = true;
            }
        }

        if ((pp.query("normalize_by_zmax", m_norm_zmax) == 1) &&
            (m_norm_zmax != 0)) {
            amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
            m_scaleFact = 1.0_rt / zmax;
            amrex::Print() << "  set scaling factor to " << m_scaleFact << '\n';
        }

    } // if forcing scheme is "indirect"
}

void ABLMesoscaleForcing::set_transition_weighting()
{
    amrex::Real zmin = m_mesh.Geom(0).ProbLo(m_axis);
    amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);

    m_blending_heights = {
        zmin, m_transition_height, m_transition_height + m_transition_thickness,
        zmax};
    m_blending_values = {1.0_rt, 1.0_rt, 0.0_rt, 0.0_rt};

    if (!m_user_specified_weighting) {
        // set weighting profile based on forcing transition
        m_weighting_heights = m_blending_heights;
        m_weighting_values = {1.0_rt, 1.0_rt, 0.0_rt, 0.0_rt};
        amrex::Print() << "setting new weighting profile" << '\n';
        for (int i = 0; i < m_weighting_heights.size(); ++i) {
            amrex::Print() << "  " << m_weighting_heights[i] << " "
                           << m_weighting_values[i] << '\n';
        }
    }
}

void ABLMesoscaleForcing::update_weights()
{
    amrex::Print() << "Updating weights" << '\n';

    for (int i = 0; i < m_nht; ++i) {
        m_W[i] =
            interp::linear(m_weighting_heights, m_weighting_values, m_zht[i]);
    }

    if (have_forcing_transition()) {
        // note the blending weights will differ from the regression weights if
        // a user-specified weighting profile is specified
        for (int i = 0; i < m_nht; ++i) {
            m_blend[i] =
                interp::linear(m_blending_heights, m_blending_values, m_zht[i]);
        }
    }
}

void ABLMesoscaleForcing::indirect_forcing_init()
{
    if (m_W.empty()) {
        // Will be here for:
        // - full profile assimilation
        // - partial profile assim w/ constant transition height
        // - partial profile assim w/ variable transition height (1st step only)
        amrex::Print() << "Initializing indirect forcing" << '\n';
        m_W.resize(m_nht);
        m_blend.resize(m_nht);
        update_weights();
    } else if (have_forcing_transition()) {
        // Will be here for:
        // - partial profile assim w/ variable transition height
        amrex::Print() << "Reinitializing indirect forcing" << '\n';
        update_weights();
    } else {
        amrex::Print() << "Should not be reinitializing indirect forcing!"
                       << '\n';
        return;
    }

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> zTz;

    // Generate the matrix Z^T W Z
    for (int irow = 0; irow < 4; irow++) {
        for (int icol = 0; icol < 4; icol++) {

            zTz(irow, icol) = 0.0_rt;

            for (int iht = 0; iht < m_nht; iht++) {
                zTz(irow, icol) =
                    zTz(irow, icol) +
                    (m_W[iht] * std::pow(
                                    m_zht[iht] * m_scaleFact,
                                    static_cast<amrex::Real>(icol + irow)));
            }
            // amrex::Print()<< "Z^T W Z ["<<irow<<","<<icol<<"] : " <<
            // zTz(irow,icol) << std::endl;
        }
    }
    // Invert the matrix Z^T W Z
    invert_mat(zTz, m_im_zTz);
}

void ABLMesoscaleForcing::invert_mat(
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
    // cppcheck-suppress constParameterReference
    amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im)
{
    amrex::Real A2323 = (m(2, 2) * m(3, 3)) - (m(2, 3) * m(3, 2));
    amrex::Real A1323 = (m(2, 1) * m(3, 3)) - (m(2, 3) * m(3, 1));
    amrex::Real A1223 = (m(2, 1) * m(3, 2)) - (m(2, 2) * m(3, 1));
    amrex::Real A0323 = (m(2, 0) * m(3, 3)) - (m(2, 3) * m(3, 0));
    amrex::Real A0223 = (m(2, 0) * m(3, 2)) - (m(2, 2) * m(3, 0));
    amrex::Real A0123 = (m(2, 0) * m(3, 1)) - (m(2, 1) * m(3, 0));
    amrex::Real A2313 = (m(1, 2) * m(3, 3)) - (m(1, 3) * m(3, 2));
    amrex::Real A1313 = (m(1, 1) * m(3, 3)) - (m(1, 3) * m(3, 1));
    amrex::Real A1213 = (m(1, 1) * m(3, 2)) - (m(1, 2) * m(3, 1));
    amrex::Real A2312 = (m(1, 2) * m(2, 3)) - (m(1, 3) * m(2, 2));
    amrex::Real A1312 = (m(1, 1) * m(2, 3)) - (m(1, 3) * m(2, 1));
    amrex::Real A1212 = (m(1, 1) * m(2, 2)) - (m(1, 2) * m(2, 1));
    amrex::Real A0313 = (m(1, 0) * m(3, 3)) - (m(1, 3) * m(3, 0));
    amrex::Real A0213 = (m(1, 0) * m(3, 2)) - (m(1, 2) * m(3, 0));
    amrex::Real A0312 = (m(1, 0) * m(2, 3)) - (m(1, 3) * m(2, 0));
    amrex::Real A0212 = (m(1, 0) * m(2, 2)) - (m(1, 2) * m(2, 0));
    amrex::Real A0113 = (m(1, 0) * m(3, 1)) - (m(1, 1) * m(3, 0));
    amrex::Real A0112 = (m(1, 0) * m(2, 1)) - (m(1, 1) * m(2, 0));

    amrex::Real det =
        (m(0, 0) * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223)) -
        (m(0, 1) * (m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223)) +
        (m(0, 2) * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123)) -
        (m(0, 3) * (m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123));
    det = 1.0_rt / det;

    im(0, 0) = det * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223);
    im(0, 1) =
        det * -((m(0, 1) * A2323) - (m(0, 2) * A1323) + (m(0, 3) * A1223));
    im(0, 2) = det * (m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213);
    im(0, 3) =
        det * -((m(0, 1) * A2312) - (m(0, 2) * A1312) + (m(0, 3) * A1212));
    im(1, 0) =
        det * -((m(1, 0) * A2323) - (m(1, 2) * A0323) + (m(1, 3) * A0223));
    im(1, 1) = det * (m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223);
    im(1, 2) =
        det * -((m(0, 0) * A2313) - (m(0, 2) * A0313) + (m(0, 3) * A0213));
    im(1, 3) = det * (m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212);
    im(2, 0) = det * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123);
    im(2, 1) =
        det * -((m(0, 0) * A1323) - (m(0, 1) * A0323) + (m(0, 3) * A0123));
    im(2, 2) = det * (m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113);
    im(2, 3) =
        det * -((m(0, 0) * A1312) - (m(0, 1) * A0312) + (m(0, 3) * A0112));
    im(3, 0) =
        det * -((m(1, 0) * A1223) - (m(1, 1) * A0223) + (m(1, 2) * A0123));
    im(3, 1) = det * (m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123);
    im(3, 2) =
        det * -((m(0, 0) * A1213) - (m(0, 1) * A0213) + (m(0, 2) * A0113));
    im(3, 3) = det * (m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112);
}

void ABLMesoscaleForcing::constant_forcing_transition(
    amrex::Vector<amrex::Real>& error)
{
    // based on SOWFA6/src/ABLForcing/drivingForce/drivingForce.C

    // find indices
    int hLevelBlend0 = -1;
    int hLevelBlend1 = -1;
    int hLevelBlendMax = -1;
    for (int iht = 0; iht < m_nht; iht++) {
        if ((hLevelBlend1 < 0) && (m_zht[iht] >= m_transition_height)) {
            hLevelBlend1 = iht;
            hLevelBlend0 = iht - 1;
        } else if (m_zht[iht] >= m_transition_height + m_transition_thickness) {
            hLevelBlendMax = iht;
            break;
        }
    }

    // error checking
    if (hLevelBlend1 < 0) {
        amrex::Print() << "Note: Did not find bottom of transition layer"
                       << '\n';
        hLevelBlend0 = m_nht - 1;
        hLevelBlend1 = m_nht - 1;
        hLevelBlendMax = m_nht - 1;
    } else if (hLevelBlendMax < 0) {
        amrex::Print() << "Note: Did not find top of transition layer" << '\n';
        hLevelBlendMax = m_nht - 1;
    }

    amrex::Print() << "Forcing transition to constant from "
                   << m_zht[hLevelBlend1] << " to " << m_zht[hLevelBlendMax]
                   << '\n';

    // calculate initial slope
    amrex::Real slope0 = (error[hLevelBlend1] - error[hLevelBlend0]) /
                         (m_zht[hLevelBlend1] - m_zht[hLevelBlend0]);
    amrex::Real dslope =
        -slope0 / (m_zht[hLevelBlendMax] - m_zht[hLevelBlend1]);

    // march from hLevelBlend1 (z >= m_transition_height)
    // up to hLevelBlendMax (z >= m_transition_height + m_transition_thickness)
    // as the slope decreases linearly to 0
    for (int iht = hLevelBlend1; iht <= hLevelBlendMax; iht++) {
        amrex::Real slope =
            slope0 + (dslope * (m_zht[iht] - m_zht[hLevelBlend1]));
        error[iht] = error[iht - 1] + (slope * (m_zht[iht] - m_zht[iht - 1]));
    }

    // set the remaining levels above hLevelBlendMax to the last value
    for (int iht = hLevelBlendMax + 1; iht < m_nht; iht++) {
        error[iht] = error[hLevelBlendMax];
    }
}

void ABLMesoscaleForcing::blend_forcings(
    const amrex::Vector<amrex::Real>& lower, // W=1
    const amrex::Vector<amrex::Real>& upper, // W=0
    amrex::Vector<amrex::Real>& error)
{
    amrex::Print() << "Blending forcings" << '\n';
    for (int iht = 0; iht < m_nht; iht++) {
        error[iht] = (m_blend[iht] * lower[iht]) +
                     ((1.0_rt - m_blend[iht]) * upper[iht]);
    }
}

} // namespace amr_wind
