#ifndef ABLOFFSHORE_TEST_UTILS_H
#define ABLOFFSHORE_TEST_UTILS_H

#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/pp_utils.H"

namespace amr_wind_tests {

// Introduce namespace for the text fixture
class ABLOffshoreMeshTest : public MeshTest
{
protected:
    void populate_parameters() override;
};

} // namespace amr_wind_tests

#endif /* ABLOFFSHORE_TEST_UTILS_H */
