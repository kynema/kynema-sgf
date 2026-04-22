#ifndef FLAT_PLATE_OPS_H
#define FLAT_PLATE_OPS_H

#include <algorithm>

#include "amr-wind/wind_energy/actuator/wing/FlatPlate.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

template <>
struct AirfoilTraits<FlatPlate>
{
    using AirfoilLookup = ThinAirfoil;
};

namespace ops {

template <>
struct ReadInputsOp<FlatPlate, ActSrcLine>
{
    void operator()(FlatPlate::DataType& data, const utils::ActParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.get_either("epsilon", wdata.eps_inp);
        pp.get("pitch", wdata.pitch);

        pp.query("chord", wdata.chord_inp);
        pp.query("cd_factor", wdata.cd_factor);

        amrex::Real max_eps = *std::ranges::max_element(wdata.eps_inp);
        amrex::Real search_radius = max_eps * wdata.chord_inp * 3.0_rt;
        const auto& p1 = wdata.start;
        const auto& p2 = wdata.end;
        // clang-format off
        info.bound_box = amrex::RealBox(
            amrex::min<amrex::Real>(p1.x(), p2.x()) - search_radius,
            amrex::min<amrex::Real>(p1.y(), p2.y()) - search_radius,
            amrex::min<amrex::Real>(p1.z(), p2.z()) - search_radius,
            amrex::max<amrex::Real>(p1.x(), p2.x()) + search_radius,
            amrex::max<amrex::Real>(p1.y(), p2.y()) + search_radius,
            amrex::max<amrex::Real>(p1.z(), p2.z()) + search_radius
        );
        // clang-format on
    }
};

template <>
struct InitDataOp<FlatPlate, ActSrcLine>
{
    void operator()(FlatPlate::DataType& data)
    {
        wing::init_data_structures(data.meta(), data.grid());

        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        wdata.chord.assign(npts, wdata.chord_inp);

        wdata.aflookup.cd_factor() = wdata.cd_factor;
    }
};

template <>
const ThinAirfoil& airfoil_lookup<FlatPlate>(FlatPlate::DataType& data)
{
    return data.meta().aflookup;
}

} // namespace ops
} // namespace amr_wind::actuator

#endif /* FLAT_PLATE_OPS_H */
