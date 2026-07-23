#ifndef DAMPINGLAYERSOURCE_H
#define DAMPINGLAYERSOURCE_H

#include "src/core/FieldRepo.H"
#include "src/equation_systems/SourceTerm.H"
#include "src/equation_systems/icns/MomentumSource.H"
#include "src/equation_systems/temperature/TemperatureSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde {

template <typename T>
class DampingLayerSource
    : public T::template Register<DampingLayerSource<T>>
{
public:
    static std::string identifier() { return "DampingLayerSource"; }

    explicit DampingLayerSource(const CFDSim& sim) { (void)sim; }

    ~DampingLayerSource() override = default;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override
    {
        (void)lev;
        (void)fstate;
        (void)src_term;
    }

private:
};

} // namespace kynema_sgf::pde

#endif /* DAMPINGLAYERSOURCE_H */
