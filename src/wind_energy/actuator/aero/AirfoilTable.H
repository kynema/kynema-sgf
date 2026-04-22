#ifndef AIRFOILTABLE_H
#define AIRFOILTABLE_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include <iosfwd>
#include <memory>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

class AirfoilLoader;

class AirfoilTable
{
public:
    friend class AirfoilLoader;

    ~AirfoilTable();

    void operator()(amrex::Real aoa, amrex::Real& cl, amrex::Real& cd) const;

    void operator()(
        amrex::Real aoa,
        amrex::Real& cl,
        amrex::Real& cd,
        amrex::Real& cm) const;

    [[nodiscard]] int num_entries() const
    {
        return static_cast<int>(m_aoa.size());
    }

    [[nodiscard]] const RealList& aoa() const { return m_aoa; }

    [[nodiscard]] const VecList& polars() const { return m_polar; }

protected:
    explicit AirfoilTable(int num_entries);

    void convert_aoa_to_radians();

    //! Angle of attack
    RealList m_aoa;

    //! Airfoil polars (Cl, Cd, Cm)
    VecList m_polar;
};

class ThinAirfoil
{
public:
    void operator()(amrex::Real aoa, amrex::Real& cl, amrex::Real& cd) const;

    amrex::Real& cd_factor() { return m_cd_factor; }

private:
    amrex::Real m_cd_factor{0.0_rt};
};

class AirfoilLoader
{
public:
    static std::unique_ptr<AirfoilTable>
    load_airfoil(const std::string& af_file, const std::string& type);

    static std::unique_ptr<AirfoilTable>
    load_text_file(const std::string& af_file);

    template <typename IStream>
    static std::unique_ptr<AirfoilTable> load_text_file(IStream& affile);

    static std::unique_ptr<AirfoilTable>
    load_openfast_airfoil(const std::string& af_file);

    template <typename IStream>
    static std::unique_ptr<AirfoilTable> load_openfast_airfoil(IStream& affile);
};

template <typename T>
struct AirfoilTraits
{
    using AirfoilLookup = AirfoilTable;
};

} // namespace amr_wind::actuator

#include "amr-wind/wind_energy/actuator/aero/AirfoilTableI.H"

#endif /* AIRFOILTABLE_H */
