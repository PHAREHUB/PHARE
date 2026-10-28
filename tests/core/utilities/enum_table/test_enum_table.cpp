#include "core/utilities/enum_table.hpp"

#include "gtest/gtest.h"

#include <cstdint>
#include <stdexcept>

using namespace PHARE::core;

namespace
{
enum class Color : std::uint16_t { red = 0, green, blue, count };

constexpr auto shuffled
    = make_enum_table<Color, int>({{Color::blue, 2}, {Color::red, 0}, {Color::green, 1}});
static_assert(shuffled[static_cast<std::size_t>(Color::red)] == 0);
static_assert(shuffled[static_cast<std::size_t>(Color::green)] == 1);
static_assert(shuffled[static_cast<std::size_t>(Color::blue)] == 2);
} // namespace

TEST(EnumTable, placesShuffledEntriesAtTheirEnumIndex)
{
    auto const table
        = make_enum_table<Color, int>({{Color::green, 1}, {Color::blue, 2}, {Color::red, 0}});

    EXPECT_EQ(table.size(), static_cast<std::size_t>(Color::count));
    EXPECT_EQ(table[static_cast<std::size_t>(Color::red)], 0);
    EXPECT_EQ(table[static_cast<std::size_t>(Color::green)], 1);
    EXPECT_EQ(table[static_cast<std::size_t>(Color::blue)], 2);
}

TEST(EnumTable, throwsOnOutOfRangeKey)
{
    EXPECT_THROW((make_enum_table<Color, int>(
                     {{Color::red, 0}, {Color::green, 1}, {Color::blue, 2}, {Color::count, 3}})),
                 std::runtime_error);
}

TEST(EnumTable, throwsOnDuplicateKey)
{
    EXPECT_THROW((make_enum_table<Color, int>(
                     {{Color::red, 0}, {Color::green, 1}, {Color::blue, 2}, {Color::red, 0}})),
                 std::runtime_error);
}

TEST(EnumTable, throwsOnMissingKey)
{
    EXPECT_THROW((make_enum_table<Color, int>({{Color::red, 0}, {Color::blue, 2}})),
                 std::runtime_error);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
