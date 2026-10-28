#include "core/models/external_field.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

#include "gtest/gtest.h"

#include <stdexcept>

using namespace PHARE::core;


TEST(ExternalField, disabledExposesNoResources)
{
    ExternalField<VecFieldMHD<2>> externalField{"ext", false};

    EXPECT_FALSE(externalField.isEnabled());
    EXPECT_TRUE(externalField.getRunTimeResourcesViewList().empty());
    EXPECT_THROW(static_cast<void>(externalField.B0()), std::runtime_error);
    EXPECT_THROW(static_cast<void>(externalField.dB0dt()), std::runtime_error);
    EXPECT_TRUE(externalField.isUsable());
}


TEST(ExternalField, enabledExposesB0AndItsTimeDerivative)
{
    ExternalField<VecFieldMHD<2>> externalField{"ext", true};
    auto& resources = externalField.getRunTimeResourcesViewList();

    EXPECT_TRUE(externalField.isEnabled());
    ASSERT_EQ(resources.size(), 2u);
    EXPECT_EQ(resources[0].name(), "ext_B0");
    EXPECT_EQ(resources[1].name(), "ext_dB0dt");
    EXPECT_EQ(&externalField.B0(), &resources[0]);
    EXPECT_EQ(&externalField.dB0dt(), &resources[1]);
    EXPECT_FALSE(externalField.isUsable());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
