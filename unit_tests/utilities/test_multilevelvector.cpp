#include "aw_test_utils/MeshTest.H"
#include "src/utilities/MultiLevelVector.H"

namespace amr_wind_tests {

class MultiLevelVectorTest : public MeshTest
{};

TEST_F(MultiLevelVectorTest, test_multilevelvector)
{
    initialize_mesh();
    amr_wind::MultiLevelVector mlv;
    mlv.resize(2, mesh().Geom());
    EXPECT_EQ(mlv.size(), 1);
    EXPECT_EQ(mlv.ncells(0), 8);
}
} // namespace amr_wind_tests
