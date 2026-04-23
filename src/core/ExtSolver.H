#ifndef EXTSOLVER_H
#define EXTSOLVER_H

#include "src/core/Factory.H"
#include "src/core/CollMgr.H"

namespace kynema_sgf {

class CFDSim;

class ExtSolver : public Factory<ExtSolver, CFDSim&>
{
public:
    static std::string base_identifier() { return "ExtSolver"; }

    ~ExtSolver() override = default;
};

class ExtSolverMgr : public CollMgr<ExtSolverMgr, ExtSolver>
{
public:
    template <typename T>
    T& get()
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<T&>(operator()(T::identifier()));
    }

    template <typename T>
    const T& get() const
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<const T&>(operator()(T::identifier()));
    }
};

} // namespace kynema_sgf

#endif /* EXTSOLVER_H */
