#include "amr-wind/incflo.H"

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("amr-wind::incflo::MakeNewLevelFromCoarse()");

    if (m_verbose > 0) {
        amrex::Print() << "Making new level " << lev << " from coarse" << '\n';
    }

    m_repo.make_new_level_from_coarse(lev, time, ba, dm);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void incflo::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("amr-wind::incflo::RemakeLevel()");

    if (m_verbose > 0) {
        amrex::Print() << "Remaking level " << lev << '\n';
    }

    m_repo.remake_level(lev, time, ba, dm);
}

// Delete level data
// overrides the pure virtual function in AmrCore
void incflo::ClearLevel(int lev)
{
    BL_PROFILE("amr-wind::incflo::ClearLevel()");
    m_repo.clear_level(lev);
}
