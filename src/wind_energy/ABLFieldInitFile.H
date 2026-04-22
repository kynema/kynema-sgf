#ifndef ABLFIELDINITFILE_H
#define ABLFIELDINITFILE_H

#include "amr-wind/core/Field.H"

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

/** Initialize subset of ABL fields using input NetCDF file
 */
class ABLFieldInitFile
{
    static_assert(AMREX_SPACEDIM == 3, "ABL requires 3 dimensional mesh");

public:
    ABLFieldInitFile();

    bool operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        int lev) const;

private:
    //! Input file with initial condition (from Machine Learning)
    std::string m_ic_input;
};

} // namespace amr_wind

#endif /* ABLFIELDINIT_H */
