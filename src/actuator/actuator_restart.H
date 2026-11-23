#ifndef ACTUATOR_RESTART_H
#define ACTUATOR_RESTART_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include "src/core/vs/vector_space.H"
#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_Utility.H"
#include "AMReX_Vector.H"

using namespace amrex::literals;

/** \file actuator_restart.H
 *
 *  Actuator restart data not held in fields (e.g., FLLC state): one text
 *  file per actuator in a directory of the CFD checkpoint. Written by the
 *  actuator root after the checkpoint fields (Actuator::write_checkpoint);
 *  read on restart from that checkpoint.
 */
namespace kynema_sgf::actuator::restart {

//! Directory of the actuators inside the checkpoint ``chkname``
inline std::string directory(const std::string& chkname)
{
    return chkname + "/actuator";
}

//! File of the actuator ``label`` inside the checkpoint ``chkname``
inline std::string
file_name(const std::string& chkname, const std::string& label)
{
    return directory(chkname) + "/" + label + ".txt";
}

/** Create the file of an actuator and write its header
 *
 *  Header: actuator label and checkpoint time. Numbers written with full
 *  round-trip precision.
 */
inline std::ofstream create_file(
    const std::string& chkname,
    const std::string& label,
    const amrex::Real time)
{
    amrex::UtilCreateDirectory(directory(chkname), 0755);
    const auto fname = file_name(chkname, label);
    std::ofstream os(fname);
    if (!os.good()) {
        amrex::Abort("Actuator " + label + ": cannot write " + fname);
    }
    os << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    os << "kynema_sgf_actuator_restart 1\n";
    os << "label " << label << '\n';
    os << "time " << time << '\n';
    return os;
}

/** Close the file of an actuator; abort when a write failed (e.g., full disk)
 *
 *  A file cut short would be refused at restart as one that does not fit.
 */
inline void close_file(
    std::ofstream& os, const std::string& chkname, const std::string& label)
{
    os.close();
    if (os.fail()) {
        amrex::Abort(
            "Actuator " + label + ": writing " + file_name(chkname, label) +
            " failed");
    }
}

/** Open the file of an actuator and read its header
 *
 *  \return False when there is no file, or it is for another actuator or
 *  time
 */
inline bool open_file(
    std::ifstream& is,
    const std::string& chkname,
    const std::string& label,
    const amrex::Real time)
{
    is.open(file_name(chkname, label));
    if (!is.good()) {
        return false;
    }
    std::string word;
    std::string label_in;
    int version = 0;
    amrex::Real time_in = 0.0_rt;
    if (!(is >> word >> version) || (word != "kynema_sgf_actuator_restart") ||
        (version != 1)) {
        return false;
    }
    if (!(is >> word >> label_in) || (word != "label") || (label_in != label)) {
        return false;
    }
    if (!(is >> word >> time_in) || (word != "time")) {
        return false;
    }
    const auto tol =
        1.0e-10_rt * amrex::max<amrex::Real>(1.0_rt, std::abs(time));
    return std::abs(time_in - time) <= tol;
}

/** Whether the checkpoint has a file of this actuator for this time
 *
 *  Written only with a state to hand over (none before the FLLC starts).
 */
inline bool has_file(
    const std::string& chkname,
    const std::string& label,
    const amrex::Real time)
{
    std::ifstream is;
    return open_file(is, chkname, label, time);
}

//! Write a named integer on a line of its own
inline void
write_value(std::ostream& os, const std::string& name, const int value)
{
    os << name << " " << value << '\n';
}

/** Read an integer written by write_value
 *
 *  \return False when the next entry has another name
 */
inline bool read_value(std::istream& is, const std::string& name, int& value)
{
    std::string name_in;
    int value_in = 0;
    if (!(is >> name_in >> value_in) || (name_in != name)) {
        return false;
    }
    value = value_in;
    return true;
}

//! Write a list of reals: name and count line, then one line per entry
template <typename RealContainer>
void write_reals(
    std::ostream& os, const std::string& name, const RealContainer& values)
{
    const int num = static_cast<int>(values.size());
    os << name << " " << num << '\n';
    for (int i = 0; i < num; ++i) {
        os << values[i] << '\n';
    }
}

/** Read a list written by write_reals
 *
 *  \param is Stream of the restart data
 *  \param name Name of the list
 *  \param values List read
 *  \param resize True to resize ``values`` to the list read, false to
 *  refuse another size
 *
 *  \return False, with ``values`` unchanged, for another name, another size
 *  (without ``resize``) or a list that ends early
 */
inline bool read_reals(
    std::istream& is,
    const std::string& name,
    amrex::Vector<amrex::Real>& values,
    const bool resize = false)
{
    std::string name_in;
    int num_in = -1;
    if (!(is >> name_in >> num_in) || (name_in != name) || (num_in < 0) ||
        (!resize && (num_in != static_cast<int>(values.size())))) {
        return false;
    }
    amrex::Vector<amrex::Real> values_in(num_in);
    for (int i = 0; i < num_in; ++i) {
        if (!(is >> values_in[i])) {
            return false;
        }
    }
    values = values_in;
    return true;
}

//! Write a list of vectors: name and count line, then one line per entry
template <typename VecContainer>
void write_list(
    std::ostream& os, const std::string& name, const VecContainer& values)
{
    const int num = static_cast<int>(values.size());
    os << name << " " << num << '\n';
    for (int i = 0; i < num; ++i) {
        os << values[i][0] << " " << values[i][1] << " " << values[i][2]
           << '\n';
    }
}

/** Read a list written by write_list into a list of the same size
 *
 *  \return False, with ``values`` unchanged, for another name, another size
 *  or a list that ends early
 */
template <typename VecContainer>
bool read_list(std::istream& is, const std::string& name, VecContainer& values)
{
    const int num = static_cast<int>(values.size());
    std::string name_in;
    int num_in = -1;
    if (!(is >> name_in >> num_in) || (name_in != name) || (num_in != num)) {
        return false;
    }
    amrex::Vector<vs::Vector> values_in(num);
    for (int i = 0; i < num; ++i) {
        if (!(is >> values_in[i][0] >> values_in[i][1] >> values_in[i][2])) {
            return false;
        }
    }
    for (int i = 0; i < num; ++i) {
        values[i] = values_in[i];
    }
    return true;
}

} // namespace kynema_sgf::actuator::restart

#endif /* ACTUATOR_RESTART_H */
