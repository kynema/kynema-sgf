#ifndef MESHTEST_H
#define MESHTEST_H

#include <memory>

#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/AmrTestMesh.H"
#include "src/core/SimTime.H"
#include "src/utilities/constants.H"
#include "src/utilities/tagging/CartBoxRefinement.H"

namespace amr_wind_tests {

/** Base class for test fixtures that require a valid AmrMesh instance
 */
class MeshTest : public AmrexTest
{
public:
    MeshTest() = default;

    ~MeshTest() override = default;

    /** Populate the ParmParse instance with all the necessary input variables
     *  for this test fixture
     */
    virtual void populate_parameters();

    //! Create the mesh instance
    virtual void create_mesh_instance();

    template <typename T>
    void create_mesh_instance()
    {
        if (!m_mesh) {
            reset_prob_domain();
            m_mesh.reset(new T());
        }
    }

    //! Perform initialization actions
    virtual void initialize_mesh();

    //! Reset the problem domain for this test
    virtual void reset_prob_domain();

    kynema_sgf::CFDSim& sim() { return mesh().sim(); }

    kynema_sgf::SimTime& time() { return sim().time(); }

    AmrTestMesh& mesh() { return *m_mesh; }

    template <typename T>
    T* mesh()
    {
        return dynamic_cast<T*>(m_mesh.get());
    }

protected:
    std::unique_ptr<AmrTestMesh> m_mesh;

    bool m_need_params{true};
};

class RefineMesh : public AmrTestMesh
{
public:
    amrex::Vector<std::unique_ptr<kynema_sgf::RefinementCriteria>>&
    refine_criteria_vec()
    {
        return m_refine_crit;
    }

protected:
    void ErrorEst(
        int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override
    {
        for (const auto& ref : m_refine_crit) {
            (*ref)(lev, tags, time, ngrow);
        }
    }

private:
    amrex::Vector<std::unique_ptr<kynema_sgf::RefinementCriteria>> m_refine_crit;
};

} // namespace amr_wind_tests

#endif /* MESHTEST_H */
