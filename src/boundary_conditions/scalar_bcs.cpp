#include "amr-wind/boundary_conditions/scalar_bcs.H"

namespace amr_wind::scalar_bc {
void register_scalar_dirichlet(
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
            "Scalar BC: Only constant dirichlet supported for Wall BC");
    }

    if (inflow_udf != "ConstDirichlet") {
        register_inflow_scalar_dirichlet<ConstDirichlet>(
            field, inflow_udf, mesh, time);
    }

    if (inflow_outflow_udf != "ConstDirichlet") {
        register_inflow_scalar_dirichlet<ConstDirichlet>(
            field, inflow_outflow_udf, mesh, time);
    }
}
} // namespace amr_wind::scalar_bc
