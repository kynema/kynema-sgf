#ifndef ABL_TEST_UTILS_H
#define ABL_TEST_UTILS_H

#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/pp_utils.H"

namespace kynema_sgf_tests {

class ABLTest : public AmrexTest
{};

// Introduce namespace for the text fixture
class ABLMeshTest : public MeshTest
{
protected:
    void populate_parameters() override;
};

} // namespace kynema_sgf_tests

#endif /* ABL_TEST_UTILS_H */
