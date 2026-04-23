#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/turbulence/TurbulenceModel.H"

namespace kynema_sgf_tests {

class TurbTest : public MeshTest
{};

TEST_F(TurbTest, test_turb_create)
{
    initialize_mesh();

    kynema_sgf::turbulence::TurbulenceModel::print(std::cout);
}

} // namespace kynema_sgf_tests
