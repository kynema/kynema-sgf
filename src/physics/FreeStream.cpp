#include "amr-wind/physics/FreeStream.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/physics/udfs/UDF.H"

namespace amr_wind {

namespace {

std::unique_ptr<udf::UDF> process_field_func(Field& fld)
{
    std::string udf_name{"ConstValue"};
    amrex::ParmParse pp("FreeStream");
    const std::string key = fld.name() + "_type";
    pp.query(key, udf_name);

    return udf::UDF::create(udf_name, fld);
}

} // namespace

FreeStream::FreeStream(const CFDSim& sim) : m_sim(sim) {}

void FreeStream::pre_init_actions()
{
    const auto& repo = m_sim.repo();
    m_field_funcs["density"] = process_field_func(repo.get_field("density"));

    const auto& pde_mgr = m_sim.pde_manager();
    {
        auto& vel = pde_mgr.icns().fields().field;
        m_field_funcs[vel.name()] = process_field_func(vel);
    }
    for (const auto& eqn : pde_mgr.scalar_eqns()) {
        auto& fld = eqn->fields().field;
        m_field_funcs[fld.name()] = process_field_func(fld);
    }

    {
        amrex::ParmParse pp("FreeStream");
        amrex::Vector<std::string> fields;
        pp.queryarr("fields", fields);
        for (const auto& fname : fields) {
            const auto it = m_field_funcs.find(fname);
            if (it != m_field_funcs.end()) {
                continue;
            }

            if (!repo.field_exists(fname)) {
                amrex::Abort("FreeStream: Invalid field requested: " + fname);
            }
            m_field_funcs[fname] = process_field_func(repo.get_field(fname));
        }
    }
}

/** Initialize the fields at the start of simulation.
 */
void FreeStream::initialize_fields(int level, const amrex::Geometry& geom)
{
    for (const auto& ff : m_field_funcs) {
        (*ff.second)(level, geom);
    }
}

} // namespace amr_wind
