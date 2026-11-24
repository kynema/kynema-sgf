#ifndef AMRTESTMESH_H
#define AMRTESTMESH_H

#include <memory>
#include <unordered_map>
#include <utility>

#include "AMReX_AmrCore.H"
#include "src/CFDSim.H"
#include "src/core/FieldRepo.H"
#include "src/core/Field.H"
#include "src/equation_systems/icns/MomentumSource.H"
#include "src/equation_systems/temperature/TemperatureSource.H"

namespace kynema_sgf_tests {

/* Base class for tests that require a valid Mesh
 *
 *  This class specializes amrex::AmrCore and provides a barebones interface
 *  that can be used within unit test fixtures.
 */
class AmrTestMesh : public amrex::AmrCore
{
public:
    AmrTestMesh();

    ~AmrTestMesh() override = default;

    // Create the initial AMR hierarchy
    virtual void initialize_mesh(amrex::Real current_time);

    // Return the total number of existing levels in AMR hierarchy
    int num_levels() const { return finest_level + 1; }

    kynema_sgf::CFDSim& sim() { return m_sim; }

    kynema_sgf::FieldRepo& field_repo() { return m_repo; }

    template <typename T, typename... Args>
    T& create_momentum_source([[maybe_unused]] Args&&... args)
    {
        auto source = std::make_unique<T>(m_sim, std::forward<Args>(args)...);
        auto& result = *source;
        m_momentum_sources.push_back(std::move(source));
        return result;
    }

    template <typename T, typename... Args>
    T& create_temperature_source([[maybe_unused]] Args&&... args)
    {
        auto source = std::make_unique<T>(m_sim, std::forward<Args>(args)...);
        auto& result = *source;
        m_temperature_sources.push_back(std::move(source));
        return result;
    }

protected:
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    void ClearLevel(int lev) override;

    void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    // Registered sources must outlive the simulation that references them.
    amrex::Vector<std::unique_ptr<kynema_sgf::pde::MomentumSource>>
        m_momentum_sources;
    amrex::Vector<std::unique_ptr<kynema_sgf::pde::TemperatureSource>>
        m_temperature_sources;
    kynema_sgf::CFDSim m_sim;
    kynema_sgf::FieldRepo& m_repo;
};

} // namespace kynema_sgf_tests

#endif /* AMRTESTMESH_H */
