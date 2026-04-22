#ifndef FREESTREAM_H
#define FREESTREAM_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/physics/udfs/UDF.H"

#include <unordered_map>
#include <memory>

namespace amr_wind {

/** Freestream physics
 *  \ingroup physics
 *
 *  Initializes the velocity, density and any specified scalar fields to
 *  user-defined constant values throughout the domain.
 */
class FreeStream : public Physics::Register<FreeStream>
{
public:
    static std::string identifier() { return "FreeStream"; }

    explicit FreeStream(const CFDSim& sim);

    ~FreeStream() override = default;

    //! Initialize the fields
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    const CFDSim& m_sim;

    std::unordered_map<std::string, std::unique_ptr<udf::UDF>> m_field_funcs;
};

} // namespace amr_wind

#endif /* FreeStream_H */
