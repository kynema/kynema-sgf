#ifndef MMS_TEST_UTILS_H
#define MMS_TEST_UTILS_H

#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/pp_utils.H"

namespace amr_wind_tests {

class MMSMeshTest : public MeshTest
{};

namespace utils {
void populate_mms_params();
}
} // namespace amr_wind_tests

#endif /* MMS_TEST_UTILS_H */
