#include "gtest/gtest.h"
#include "src/equation_systems/source_terms/DampingLayerSource.H"

namespace kynema_sgf_tests {

TEST(DampingLayerSource, target_type_parse)
{
    EXPECT_EQ(string_to_target_type("constant"), TargetType::Constant);
    EXPECT_EQ(string_to_target_type("PROFILE"), TargetType::Profile);
    EXPECT_EQ(string_to_target_type("Function"), TargetType::Function);
    EXPECT_EQ(string_to_target_type("field"), TargetType::Field);
}

TEST(DampingLayerSource, trait_field_names)
{
    using namespace kynema_sgf::pde;

    EXPECT_STREQ(
        DampingLayerSourceTraits<MomentumSource>::field_name, "velocity");
    EXPECT_STREQ(
        DampingLayerSourceTraits<TemperatureSource>::field_name, "temperature");
    EXPECT_STREQ(
        DampingLayerSourceTraits<DensitySource>::field_name, "density");
    EXPECT_STREQ(DampingLayerSourceTraits<TKESource>::field_name, "tke");
    EXPECT_STREQ(DampingLayerSourceTraits<SDRSource>::field_name, "sdr");
}

} // namespace kynema_sgf_tests
