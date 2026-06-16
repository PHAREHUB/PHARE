#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "core/data/vector.hpp"

using namespace PHARE::core;


TEST(MinimizingVector, defaultConstruct)
{
    MinimizingVector<int> v;
    EXPECT_EQ(v.size(), 0u);
    EXPECT_EQ(v.capacity(), 0u);
}


TEST(MinimizingVector, getResizesVector)
{
    MinimizingVector<int> v;
    v.get(10);
    EXPECT_EQ(v.size(), 10u);
}


TEST(MinimizingVector, getPreservesContent)
{
    MinimizingVector<int> v;
    v.get(5);
    for (std::size_t i = 0; i < 5; ++i)
        v[i] = static_cast<int>(i);

    v.get(5);
    for (std::size_t i = 0; i < 5; ++i)
        EXPECT_EQ(v[i], static_cast<int>(i));
}


TEST(MinimizingVector, getNoNoCopyResizesVector)
{
    MinimizingVector<int> v;
    v.get_no_copy(20);
    EXPECT_EQ(v.size(), 20u);
}


TEST(MinimizingVector, reserveAndClearSetsCapacityAndZeroSize)
{
    MinimizingVector<int> v;
    v.reserve_and_clear(50);
    EXPECT_EQ(v.size(), 0u);
    EXPECT_GE(v.capacity(), 50u);
}


TEST(MinimizingVector, clearSetsSize0)
{
    MinimizingVector<int> v;
    v.get(10);
    v.clear();
    EXPECT_EQ(v.size(), 0u);
}


TEST(MinimizingVector, destroyReleasesMemory)
{
    MinimizingVector<int> v;
    v.get(100);
    v.destroy();
    EXPECT_EQ(v.size(), 0u);
    EXPECT_EQ(v.capacity(), 0u);
}


TEST(MinimizingVector, resizeIsAliasForGet)
{
    MinimizingVector<int> v;
    v.resize(7);
    EXPECT_EQ(v.size(), 7u);
}


TEST(MinimizingVector, emplaceBack)
{
    MinimizingVector<int> v;
    v.emplace_back(42);
    v.emplace_back(99);
    EXPECT_EQ(v.size(), 2u);
    EXPECT_EQ(v[0], 42);
    EXPECT_EQ(v[1], 99);
}


TEST(MinimizingVector, pushBack)
{
    MinimizingVector<int> v;
    v.push_back(7);
    v.push_back(13);
    EXPECT_EQ(v.size(), 2u);
    EXPECT_EQ(v[0], 7);
    EXPECT_EQ(v[1], 13);
}


TEST(MinimizingVector, dataPointerMatchesElements)
{
    MinimizingVector<int> v;
    v.get(4);
    v[0] = 1; v[1] = 2; v[2] = 3; v[3] = 4;
    auto* p = v.data();
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(p[i], i + 1);
}


TEST(MinimizingVector, operatorStarReturnsUnderlyingVector)
{
    MinimizingVector<int> v;
    v.get(3);
    v[0] = 10; v[1] = 20; v[2] = 30;
    auto& inner = *v;
    EXPECT_EQ(inner.size(), 3u);
    EXPECT_EQ(inner[1], 20);
}


TEST(MinimizingVector, operatorCallReturnsUnderlyingVector)
{
    MinimizingVector<int> v;
    v.get(3);
    v[0] = 5;
    EXPECT_EQ(v()[0], 5);
}


// After `period` consecutive accesses where s < capacity * percentile,
// the vector should shrink its capacity toward capacity * realloc_to.
TEST(MinimizingVector, shrinksAfterPeriodSmallRequests)
{
    MinimizingVector<int> v;
    // Establish a large capacity
    v.get(1000);
    auto const large_cap = v.capacity();
    EXPECT_GE(large_cap, 1000u);

    // Drive period (100) consecutive small requests — each well below percentile (80%)
    std::size_t const small = 100; // 100 << 800 (= 1000 * 0.8)
    for (std::size_t i = 0; i < v.period; ++i)
        v.get(small);

    // After exactly period requests the shrink fires
    EXPECT_LT(v.capacity(), large_cap);
}


TEST(MinimizingVector, doesNotShrinkBeforePeriod)
{
    MinimizingVector<int> v;
    v.get(1000);
    auto const large_cap = v.capacity();

    std::size_t const small = 100;
    for (std::size_t i = 0; i < v.period - 1; ++i)
        v.get(small);

    EXPECT_EQ(v.capacity(), large_cap);
}


TEST(MinimizingVector, counterResetWhenLargeRequest)
{
    MinimizingVector<int> v;
    v.get(1000);
    auto const large_cap = v.capacity();

    // 50 small requests, then one large — counter should reset
    for (std::size_t i = 0; i < 50; ++i)
        v.get(100);

    v.get(1000); // large — resets counter

    // 50 more small requests: total small streak is 50, not 100 → no shrink yet
    for (std::size_t i = 0; i < 50; ++i)
        v.get(100);

    EXPECT_EQ(v.capacity(), large_cap);
}


TEST(MinimizingVector, getNoCopyShrinks)
{
    MinimizingVector<int> v;
    v.get_no_copy(1000);
    auto const large_cap = v.capacity();

    for (std::size_t i = 0; i < v.period; ++i)
        v.get_no_copy(100);

    EXPECT_LT(v.capacity(), large_cap);
}


TEST(MinimizingVector, reserveAndClearShrinks)
{
    MinimizingVector<int> v;
    v.reserve_and_clear(1000);
    auto const large_cap = v.capacity();

    for (std::size_t i = 0; i < v.period; ++i)
        v.reserve_and_clear(100);

    EXPECT_LT(v.capacity(), large_cap);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
