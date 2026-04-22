#ifndef VELOCITY_BCS_H
#define VELOCITY_BCS_H

#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/physics/udfs/LinearProfile.H"
#include "amr-wind/physics/udfs/PowerLawProfile.H"
#include "amr-wind/physics/udfs/BurggrafLid.H"
#include "amr-wind/physics/udfs/Rankine.H"
#include "amr-wind/physics/udfs/CustomVelocity.H"
#include "amr-wind/physics/udfs/TwoLayer.H"

namespace amr_wind::vel_bc {

template <typename WallOp>
void register_inflow_vel_dirichlet(
    Field& field,
    const std::string& inflow_udf,
    const amrex::AmrCore& mesh,
    const SimTime& time)
{
    if (inflow_udf == "LinearProfile") {
        using InflowOp = BCOpCreator<udf::LinearProfile, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else if (inflow_udf == "PowerLawProfile") {
        using InflowOp = BCOpCreator<udf::PowerLawProfile, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else if (inflow_udf == "BurggrafLid") {
        using InflowOp = BCOpCreator<udf::BurggrafLid, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else if (inflow_udf == "Rankine") {
        using InflowOp = BCOpCreator<udf::Rankine, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else if (inflow_udf == "CustomVelocity") {
        using InflowOp = BCOpCreator<udf::CustomVelocity, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else if (inflow_udf == "TwoLayer") {
        using InflowOp = BCOpCreator<udf::TwoLayer, WallOp>;
        field.register_fill_patch_op<FieldFillPatchOps<InflowOp>>(
            mesh, time, InflowOp(field));
    } else {
        amrex::Abort("Velocity BC: Invalid dirichlet BC type = " + inflow_udf);
    }
}

void register_velocity_dirichlet(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    const amrex::Array<const std::string, 3>& udfs)
{
    const std::string& inflow_udf = udfs[0];
    const std::string& inflow_outflow_udf = udfs[1];
    const std::string& wall_udf = udfs[2];

    if ((inflow_udf == "ConstDirichlet") &&
        (inflow_outflow_udf == "ConstDirichlet") &&
        (wall_udf == "ConstDirichlet")) {
        return;
    }

    if (wall_udf != "ConstDirichlet") {
        amrex::Abort(
            "Velocity BC: Only constant dirichlet supported for Wall BC");
    }

    if (inflow_udf != "ConstDirichlet") {
        register_inflow_vel_dirichlet<ConstDirichlet>(
            field, inflow_udf, mesh, time);
    }

    if (inflow_outflow_udf != "ConstDirichlet") {
        register_inflow_vel_dirichlet<ConstDirichlet>(
            field, inflow_outflow_udf, mesh, time);
    }
}

} // namespace amr_wind::vel_bc

#endif /* VELOCITY_BCS_H */
