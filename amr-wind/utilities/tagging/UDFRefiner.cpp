#include "amr-wind/utilities/tagging/UDFRefiner.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::tagging {

UDFRefiner::UDFRefiner(const CFDSim& sim, const std::string& key)
    : m_sim(sim), m_bound_box(m_sim.repo().mesh().Geom(0).ProbDomain())
{
    amrex::ParmParse pp(key);
    std::string udf_str;
    pp.get("udf", udf_str);

    amrex::Vector<amrex::Real> box_lo(AMREX_SPACEDIM, 0);
    amrex::Vector<amrex::Real> box_hi(AMREX_SPACEDIM, 0);
    if (pp.queryarr("box_lo", box_lo, 0, static_cast<int>(box_lo.size())) ==
        1) {
        m_bound_box.setLo(box_lo);
    }
    if (pp.queryarr("box_hi", box_hi, 0, static_cast<int>(box_hi.size())) ==
        1) {
        m_bound_box.setHi(box_hi);
    }
    m_parser.define(udf_str);
    m_parser.registerVariables({"t", "x", "y", "z"});
}

void UDFRefiner::operator()(
    const amrex::Box& bx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::TagBox::TagType>& tag) const
{
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const auto time = m_sim.time().current_time();
    auto udf_func = m_parser.compile<AMREX_SPACEDIM + 1>();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
        const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
        const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
        const auto val = udf_func(time, x, y, z);
        if (static_cast<bool>(val)) {
            tag(i, j, k) = amrex::TagBox::SET;
        }
    });
}

} // namespace amr_wind::tagging
