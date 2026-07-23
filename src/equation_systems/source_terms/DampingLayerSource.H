#ifndef DAMPINGLAYERSOURCE_H
#define DAMPINGLAYERSOURCE_H

#include <string>

#include "src/core/FieldRepo.H"
#include "src/equation_systems/density/DensitySource.H"
#include "src/equation_systems/sdr/SDRSource.H"
#include "src/equation_systems/SourceTerm.H"
#include "src/equation_systems/icns/MomentumSource.H"
#include "src/equation_systems/temperature/TemperatureSource.H"
#include "src/equation_systems/tke/TKESource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde {

template <typename T>
struct DampingLayerSourceTraits
{
    static constexpr bool is_specialized = false;
    static constexpr const char* field_name = "";
};

template <>
struct DampingLayerSourceTraits<MomentumSource>
{
    static constexpr bool is_specialized = true;
    static constexpr const char* field_name = "velocity";
};

template <>
struct DampingLayerSourceTraits<TemperatureSource>
{
    static constexpr bool is_specialized = true;
    static constexpr const char* field_name = "temperature";
};

template <>
struct DampingLayerSourceTraits<DensitySource>
{
    static constexpr bool is_specialized = true;
    static constexpr const char* field_name = "density";
};

template <>
struct DampingLayerSourceTraits<TKESource>
{
    static constexpr bool is_specialized = true;
    static constexpr const char* field_name = "tke";
};

template <>
struct DampingLayerSourceTraits<SDRSource>
{
    static constexpr bool is_specialized = true;
    static constexpr const char* field_name = "sdr";
};

template <>
struct DampingLayerSourceTraits<SourceTerm>
{
    static constexpr bool is_specialized = true;
    static constexpr const char* field_name = "";
};

template <typename T>
class DampingLayerSource
    : public T::template Register<DampingLayerSource<T>>
{
public:
    static_assert(
        DampingLayerSourceTraits<T>::is_specialized,
        "DampingLayerSourceTraits<T> specialization missing for this source "
        "type. Add a specialization with field_name.");

    static std::string identifier()
    {
        return "DampingLayerSource";
    }

    static constexpr const char* target_field_name()
    {
        return DampingLayerSourceTraits<T>::field_name;
    }

    explicit DampingLayerSource(const CFDSim& sim)
        : m_target_field_name(target_field_name())
    {
        (void)sim;
    }

    ~DampingLayerSource() override = default;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override
    {
        (void)lev;
        (void)fstate;
        (void)src_term;
    }

private:
    std::string m_target_field_name;
};

} // namespace kynema_sgf::pde

#endif /* DAMPINGLAYERSOURCE_H */
