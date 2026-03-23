#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include <fstream>
#include <algorithm>

namespace amr_wind::actuator {

AirfoilTable::AirfoilTable(const int num_entries)
    : m_aoa(num_entries), m_polar(num_entries)
{}

AirfoilTable::~AirfoilTable() = default;

void AirfoilTable::operator()(
    const amrex::Real aoa, amrex::Real& cl, amrex::Real& cd) const
{
    namespace interp = ::amr_wind::interp;
    vs::Vector polar = interp::linear(m_aoa, m_polar, aoa);
    cl = polar.x();
    cd = polar.y();
}

void AirfoilTable::operator()(
    const amrex::Real aoa,
    amrex::Real& cl,
    amrex::Real& cd,
    amrex::Real& cm) const
{
    namespace interp = ::amr_wind::interp;
    vs::Vector polar = interp::linear(m_aoa, m_polar, aoa);
    cl = polar.x();
    cd = polar.y();
    cm = polar.z();
}

void ThinAirfoil::operator()(
    const amrex::Real aoa, amrex::Real& cl, amrex::Real& cd) const
{
    cl = ::amr_wind::utils::two_pi() * aoa;
    cd = m_cd_factor * std::sin(aoa);
}

void AirfoilTable::convert_aoa_to_radians()
{
    std::ranges::transform(m_aoa, m_aoa.begin(), [](amrex::Real aoa_in) {
        return utils::radians(aoa_in);
    });
}

std::unique_ptr<AirfoilTable>
AirfoilLoader::load_text_file(const std::string& af_file)
{
    std::ifstream afh(af_file);

    if (!afh.good()) {
        amrex::Abort("AirfoilLoader: Cannot open airfoil file: " + af_file);
    }
    return load_text_file(afh);
}

std::unique_ptr<AirfoilTable>
AirfoilLoader::load_openfast_airfoil(const std::string& af_file)
{
    std::ifstream afh(af_file);

    if (!afh.good()) {
        amrex::Abort("AirfoilLoader: Cannot open airfoil file: " + af_file);
    }
    return load_openfast_airfoil(afh);
}

std::unique_ptr<AirfoilTable>
AirfoilLoader::load_airfoil(const std::string& af_file, const std::string& type)
{
    const std::string& aftype = amrex::toLower(type);

    std::unique_ptr<AirfoilTable> af;
    if (aftype == "openfast") {
        af = load_openfast_airfoil(af_file);
    } else if (aftype == "text") {
        af = load_text_file(af_file);
    } else {
        amrex::Abort("Invalid airfoil type specified");
    }

    return af;
}

} // namespace amr_wind::actuator
