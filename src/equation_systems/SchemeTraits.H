#ifndef SCHEMETRAITS_H
#define SCHEMETRAITS_H

#include <string>

namespace amr_wind::fvm {

/** Characteristics of the Godnov numerical scheme
 *  \ingroup eqsys
 */
struct Godunov
{
    //! String used to uniquely identify this scheme
    static std::string scheme_name() { return "Godunov"; }

    //! Number of states in the field variable
    static constexpr int num_states = 2;
    //! Number of ghost in the state variable
    static constexpr int nghost_state = 3;
    //! Number of ghost cells in the source term variable
    static constexpr int nghost_src = 1;
    //! Number of ghost cells in the MAC face variables
    static constexpr int nghost_mac = 1;

    //! Number of states for diffusion field
    static constexpr int num_diff_states = 1;
    //! Number of states for convective term
    static constexpr int num_conv_states = 1;
};

/** Characteristics of the MOL numerical scheme
 *  \ingroup eqsys
 */
struct MOL
{
    static std::string scheme_name() { return "MOL"; }

    //! Number of states in the field variable
    static constexpr int num_states = 2;
    //! Number of ghost cells in the  state variable
    static constexpr int nghost_state = 2;
    //! Number of ghost cells in the source term variable
    static constexpr int nghost_src = 0;
    //! Number of ghost cells in the MAC face variables
    static constexpr int nghost_mac = 0;

    //! Number of states for diffusion field
    static constexpr int num_diff_states = 2;
    //! Number of states for convective term
    static constexpr int num_conv_states = 2;
};

} // namespace amr_wind::fvm

#endif /* SCHEMETRAITS_H */
