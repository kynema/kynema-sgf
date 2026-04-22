#include "AmrTestMesh.H"
#include "gtest/gtest.h"
#include "AMReX_MultiFab.H"

namespace amr_wind_tests {

AmrTestMesh::AmrTestMesh() : m_sim(*this), m_repo(m_sim.repo())
{
    m_sim.time().parse_parameters();
}

void AmrTestMesh::initialize_mesh(amrex::Real current_time)
{
    InitFromScratch(current_time);
}

void AmrTestMesh::MakeNewLevelFromScratch(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    m_repo.make_new_level_from_scratch(lev, time, ba, dm);
}

void AmrTestMesh::MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // Concrete test fixtures must override and initialize field data if
    // necessary

    m_repo.make_new_level_from_coarse(lev, time, ba, dm);
}

void AmrTestMesh::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    m_repo.remake_level(lev, time, ba, dm);
}

void AmrTestMesh::ClearLevel(int lev) { m_repo.clear_level(lev); }

void AmrTestMesh::ErrorEst(
    int /* lev */,
    amrex::TagBoxArray& /* tags */,
    amrex::Real /* time */,
    int /* ngrow */)
{
    amrex::Abort("Not implemented");
}

} // namespace amr_wind_tests
