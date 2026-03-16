#ifndef AIRFOILTABLEI_H
#define AIRFOILTABLEI_H

#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"

#include <string>

namespace amr_wind::actuator {

template <typename IStream>
std::unique_ptr<AirfoilTable> AirfoilLoader::load_text_file(IStream& affile)
{
    int num_entries;
    affile >> num_entries;

    std::unique_ptr<AirfoilTable> aftab(new AirfoilTable(num_entries));
    for (int i = 0; i < num_entries; ++i) {
        auto& pp = aftab->m_polar[i];
        affile >> aftab->m_aoa[i] >> pp.x() >> pp.y() >> pp.z();
    }

    aftab->convert_aoa_to_radians();
    return aftab;
}

template <typename IStream>
std::unique_ptr<AirfoilTable>
AirfoilLoader::load_openfast_airfoil(IStream& affile)
{
    std::string buf;
    int num_entries = -1;
    while (affile.good() && num_entries < 0) {
        std::getline(affile, buf);
        const auto found = buf.find("NumAlf");
        if (found != std::string::npos) {
            std::stringstream ss(buf);
            ss >> num_entries;
        }
    }

    if (!affile.good() && (num_entries < 0)) {
        amrex::Abort("AirfoilLoader: Error reading OpenFAST airfoil file");
    }

    // Skip two comment lines
    std::getline(affile, buf);
    std::getline(affile, buf);

    std::unique_ptr<AirfoilTable> aftab(new AirfoilTable(num_entries));
    for (int i = 0; i < num_entries; ++i) {
        auto& pp = aftab->m_polar[i];
        affile >> aftab->m_aoa[i] >> pp.x() >> pp.y() >> pp.z();
    }

    aftab->convert_aoa_to_radians();
    return aftab;
}

} // namespace amr_wind::actuator
#endif
