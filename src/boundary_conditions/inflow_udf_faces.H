#ifndef INFLOW_UDF_FACES_H
#define INFLOW_UDF_FACES_H

#include "src/core/Field.H"
#include "src/incflo_enums.H"

#include "AMReX_ParmParse.H"

#include <string>

namespace kynema_sgf::bc_udf {

/** Check that the one inflow UDF of a field can fill all its inflow faces
 *
 *  A single fill operator serves every mass_inflow and mass_inflow_outflow
 *  face of a field, and it calls the inflow UDF on each of them. So the two
 *  kinds of face cannot use different UDFs, and a UDF selected on some faces
 *  would silently replace the constant value of the others. TabulatedProfile
 *  is the exception: it reads the key of each face and gives a face that does
 *  not select it its constant value.
 *
 *  \param field Field whose inflow faces are checked
 *  \param udfs UDFs of its mass_inflow, mass_inflow_outflow and wall faces
 *  \param label Start of the error message, such as "Velocity BC"
 */
inline void check_inflow_udf_faces(
    const Field& field,
    const amrex::Array<const std::string, 3>& udfs,
    const std::string& label)
{
    const std::string& inflow_udf = udfs[0];
    const std::string& inflow_outflow_udf = udfs[1];
    const std::string& fname = field.name();

    if ((inflow_udf != "ConstDirichlet") &&
        (inflow_outflow_udf != "ConstDirichlet") &&
        (inflow_udf != inflow_outflow_udf)) {
        amrex::Abort(
            label + ": " + fname + ".inflow_type = " + inflow_udf +
            " (mass_inflow faces) and " + fname +
            ".inflow_outflow_type = " + inflow_outflow_udf +
            " (mass_inflow_outflow faces) differ; one UDF fills every "
            "inflow face, so use the same one on both");
    }

    const std::string& udf =
        (inflow_udf != "ConstDirichlet") ? inflow_udf : inflow_outflow_udf;
    if ((udf == "ConstDirichlet") || (udf == "TabulatedProfile")) {
        return;
    }

    const amrex::Array<std::string, 6> face_names{"xlo", "ylo", "zlo",
                                                  "xhi", "yhi", "zhi"};
    const auto& bctype = field.bc_type();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        const auto ori = oit();
        std::string key;
        if (bctype[ori] == BC::mass_inflow) {
            key = fname + ".inflow_type";
        } else if (bctype[ori] == BC::mass_inflow_outflow) {
            key = fname + ".inflow_outflow_type";
        } else {
            continue;
        }
        std::string selected{"ConstDirichlet"};
        amrex::ParmParse(face_names[ori]).query(key, selected);
        if (selected != udf) {
            amrex::Abort(
                label + ": " + face_names[ori] + "." + key + " is " + selected +
                ", but " + udf + " fills every inflow face of " + fname +
                " (only TabulatedProfile falls back to the constant face by "
                "face); select " +
                udf + " on this face too");
        }
    }
}

} // namespace kynema_sgf::bc_udf

#endif /* INFLOW_UDF_FACES_H */
