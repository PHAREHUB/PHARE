
#include <string>

#include "utilities/range/range.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"


using PHARE::makeRange;
using PHARE::Range;

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




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
