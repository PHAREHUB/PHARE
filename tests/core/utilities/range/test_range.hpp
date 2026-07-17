#ifndef PHARE_TEST_UTILITY_TEST_RANGE_HPP
#define PHARE_TEST_UTILITY_TEST_RANGE_HPP

#include <string>

#include "core/utilities/range/range.hpp"


#include "gmock/gmock.h"
#include "gtest/gtest.h"


using PHARE::core::makeRange;
using PHARE::core::Range;

TEST(ARangeOnEmptyContainer, hasZeroSize)
{
    std::vector<double> b;
    auto range = makeRange(std::begin(b), std::end(b));
    EXPECT_EQ(0, range.size());
}


TEST(ARange, canReturnSize)
{
    std::vector<double> b(10);
    auto range = makeRange(std::begin(b), std::end(b));
    EXPECT_EQ(10, range.size());
}



TEST(ARange, returnBeginAndEnd)
{
    std::vector<double> b(10);
    auto range = makeRange(std::begin(b), std::end(b));
    EXPECT_EQ(std::begin(b), range.begin());
    EXPECT_EQ(std::end(b), range.end());
}


#endif /*  PHARE_TEST_UTILITY_TEST_RANGE_HPP  */
