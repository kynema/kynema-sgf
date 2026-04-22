#ifndef CYLINDER_H
#define CYLINDER_H

#include "src/immersed_boundary/bluff_body/BluffBody.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::ib {

struct CylinderData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real radius{0.1_rt};
    amrex::Real height{1.0_rt};
};

struct Cylinder : public BluffBodyType
{
    using InfoType = IBInfo;
    using MetaType = CylinderData;
    using DataType = IBDataHolder<Cylinder>;

    static std::string identifier() { return "Cylinder"; }
};

} // namespace kynema_sgf::ib

#endif /* CYLINDER_H */
