
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_range.hpp"
#include "test_ranges.hpp"
#include "test_range_replacer.hpp"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
