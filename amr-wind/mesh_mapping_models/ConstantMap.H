#ifndef CONSTANTMAP_H
#define CONSTANTMAP_H

#include "amr-wind/core/MeshMap.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::const_map {

/** Coonstant scaling mesh map
 *  \ingroup mesh_map
 */
class ConstantMap : public MeshMap::Register<ConstantMap>
{
public:
    static std::string identifier() { return "ConstantMap"; }

    explicit ConstantMap();

    ~ConstantMap() override = default;

    //! Construct the mesh scaling field
    void create_map(int /*lev*/, const amrex::Geometry& /*geom*/) override;

    //! Construct mesh scaling field on cell centers and nodes
    void create_cell_node_map(int /*lev*/);

    //! Construct mesh scaling field on cell faces
    void create_face_map(int /*lev*/);

    //! Construct the non-uniform mesh field
    void create_non_uniform_mesh(int /*lev*/, const amrex::Geometry& /*geom*/);

private:
    //! Factor to scale the mesh by
    amrex::Vector<amrex::Real> m_fac{1.0_rt, 1.0_rt, 1.0_rt};
};

} // namespace amr_wind::const_map

#endif /* CONSTANTMAP_H */
