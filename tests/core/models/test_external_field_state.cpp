#include "core/data/field/zero_field.hpp"
#include "core/models/external_field_state.hpp"

#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

#include "gtest/gtest.h"

#include <stdexcept>
#include <type_traits>

using namespace PHARE::core;

namespace
{
template<typename T>
constexpr bool is_zero_v = std::is_same_v<std::decay_t<T>, ZeroVecField<VecFieldMHD<2>>>;

auto const isZero = [](auto const& vecfield) { return is_zero_v<decltype(vecfield)>; };

auto isTheResource(VecFieldMHD<2> const& resource)
{
    return [&resource](auto const& vecfield) {
        if constexpr (is_zero_v<decltype(vecfield)>)
            return false;
        else
            return &vecfield == &resource;
    };
}
} // namespace


TEST(ExternalFieldState, disabledExposesNoResources)
{
    ExternalFieldState<VecFieldMHD<2>> externalField{"ext", false};

    EXPECT_FALSE(externalField.isEnabled());
    EXPECT_FALSE(externalField.isTimeDependent());
    EXPECT_TRUE(externalField.getRunTimeResourcesViewList().empty());
    EXPECT_THROW(static_cast<void>(externalField.B0()), std::runtime_error);
    EXPECT_THROW(static_cast<void>(externalField.dB0dt()), std::runtime_error);
    EXPECT_TRUE(externalField.visitB0(isZero));
    EXPECT_TRUE(externalField.visitDB0dt(isZero));
    EXPECT_TRUE(externalField.isUsable());
}


TEST(ExternalFieldState, staticExposesB0Only)
{
    ExternalFieldState<VecFieldMHD<2>> externalField{"ext", true};
    auto& resources = externalField.getRunTimeResourcesViewList();

    EXPECT_TRUE(externalField.isEnabled());
    EXPECT_FALSE(externalField.isTimeDependent());
    ASSERT_EQ(resources.size(), 1u);
    EXPECT_EQ(resources[0].name(), "ext_B0");
    EXPECT_EQ(&externalField.B0(), &resources[0]);
    EXPECT_THROW(static_cast<void>(externalField.dB0dt()), std::runtime_error);
    EXPECT_TRUE(externalField.visitB0(isTheResource(resources[0])));
    EXPECT_TRUE(externalField.visitDB0dt(isZero));
    EXPECT_FALSE(externalField.isUsable());
}


TEST(ExternalFieldState, timeDependentExposesB0AndItsTimeDerivative)
{
    ExternalFieldState<VecFieldMHD<2>> externalField{"ext", true, true};
    auto& resources = externalField.getRunTimeResourcesViewList();

    EXPECT_TRUE(externalField.isEnabled());
    EXPECT_TRUE(externalField.isTimeDependent());
    ASSERT_EQ(resources.size(), 2u);
    EXPECT_EQ(resources[0].name(), "ext_B0");
    EXPECT_EQ(resources[1].name(), "ext_dB0dt");
    EXPECT_EQ(&externalField.B0(), &resources[0]);
    EXPECT_EQ(&externalField.dB0dt(), &resources[1]);
    EXPECT_TRUE(externalField.visitB0(isTheResource(resources[0])));
    EXPECT_TRUE(externalField.visitDB0dt(isTheResource(resources[1])));
    EXPECT_FALSE(externalField.isUsable());
}


TEST(ExternalFieldState, timeDependentWithoutFieldThrows)
{
    EXPECT_THROW((ExternalFieldState<VecFieldMHD<2>>{"ext", false, true}), std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
