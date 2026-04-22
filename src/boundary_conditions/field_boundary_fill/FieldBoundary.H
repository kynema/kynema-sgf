#ifndef FIELDBOUNDARY_H
#define FIELDBOUNDARY_H

#include "src/core/Factory.H"
#include "src/core/CollMgr.H"
#include "src/core/Field.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"

namespace kynema_sgf {

class CFDSim;

/** Abstract representation of field-boundary fill models
 *
 *  This class mirrors the runtime-registration model used by Physics and
 *  provides hook points for boundary-fill specific setup and time-stepping
 *  actions.
 */
class FieldBoundary : public Factory<FieldBoundary, CFDSim&>
{
public:
    static std::string base_identifier() { return "FieldBoundary"; }

    ~FieldBoundary() override = default;

    //! Actions to perform post initialization
    virtual void post_init_actions() = 0;

    //! Perform tasks necessary before advancing timestep
    virtual void pre_advance_work() = 0;

    //! Perform tasks necessary only once per timestep, after pre_advance
    virtual void pre_predictor_work() {}

    //! Perform tasks necessary after advancing timestep
    virtual void post_advance_work() = 0;

    //! Velocity fill function to be used by nodal projection
    virtual void set_velocity(
        int /*lev*/,
        amrex::Real /*time*/,
        const Field& /*vel_fld*/,
        amrex::MultiFab& /*vel_mfab*/,
        int /*dcomp*/ = 0,
        int /*orig_comp*/ = 0) const = 0;
};

/** A collection of FieldBoundary instances active during a simulation */
class FieldBoundaryMgr : public CollMgr<FieldBoundaryMgr, FieldBoundary>
{
public:
    //! Return a concrete field-boundary instance
    template <typename T>
    T& get()
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<T&>(operator()(T::identifier()));
    }

    //! Return a concrete field-boundary instance
    template <typename T>
    [[nodiscard]] const T& get() const
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<const T&>(operator()(T::identifier()));
    }
};

} // namespace kynema_sgf

#endif /* FIELDBOUNDARY_H */
