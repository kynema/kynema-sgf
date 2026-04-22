#ifndef SCALAR_BCS_H
#define SCALAR_BCS_H

#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/physics/udfs/CustomScalar.H"

namespace amr_wind::scalar_bc {

template <typename WallOp>
void register_inflow_scalar_dirichlet(
    Field& field,
    const std::string& inflow_udf,
    const amrex::AmrCore& mesh,
    const SimTime& time)
{
    if (inflow_udf == "CustomScalar") {
        using InflowOp = BCOpCreator<udf::CustomScalar, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else {
        amrex::Abort("Scalar BC: Invalid dirichlet BC type = " + inflow_udf);
    }
}

void register_scalar_dirichlet(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    const amrex::Array<const std::string, 3>& udfs);

} // namespace amr_wind::scalar_bc

#endif /* SCALAR_BCS_H */
