#ifndef SPHERE_H
#define SPHERE_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ib {

struct SphereData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real radius{0.1_rt};
};

struct Sphere : public BluffBodyType
{
    using InfoType = IBInfo;
    using MetaType = SphereData;
    using DataType = IBDataHolder<Sphere>;

    static std::string identifier() { return "Sphere"; }
};

} // namespace amr_wind::ib

#endif /* Sphere_H */
