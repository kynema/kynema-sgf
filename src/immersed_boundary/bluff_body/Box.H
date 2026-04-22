#ifndef BOX_H
#define BOX_H

#include "src/immersed_boundary/bluff_body/BluffBody.H"

namespace kynema_sgf::ib {

struct BoxData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real length{0.1_rt};
    amrex::Real width{0.1_rt};
    amrex::Real height{0.1_rt};
};

struct Box : public BluffBodyType
{
    using InfoType = IBInfo;
    using MetaType = BoxData;
    using DataType = IBDataHolder<Box>;

    static std::string identifier() { return "Box"; }
};

} // namespace kynema_sgf::ib

#endif /* BOX_H */
