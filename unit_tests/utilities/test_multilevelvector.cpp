#include "ks_test_utils/MeshTest.H"
#include "src/utilities/MultiLevelVector.H"

namespace kynema_sgf_tests {

class MultiLevelVectorTest : public MeshTest
{};

TEST_F(MultiLevelVectorTest, test_multilevelvector)
{
    initialize_mesh();
    kynema_sgf::MultiLevelVector mlv;
    mlv.resize(2, mesh().Geom());
    EXPECT_EQ(mlv.size(), 1);
    EXPECT_EQ(mlv.ncells(0), 8);
}
} // namespace kynema_sgf_tests
