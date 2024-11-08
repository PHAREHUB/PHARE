#include <ctype.h>
#include <string>

#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"
#include "core/data/tiles/tiles.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>


using namespace PHARE::core;


template<typename TileSet>
class TileTest : public ::testing::Test
{
public:
    TileTest() {}
};

template<std::size_t dim>
class TileMock
{
public:
    static auto constexpr dimension = dim;
    Point<int, dim> lower;
    Point<int, dim> upper;


    auto size() const { return Box<int, dim>{lower, upper}.size(); }
};

using DimTiles = testing::Types<TileSet<TileMock<1>>, TileSet<TileMock<2>>, TileSet<TileMock<3>>>;

TYPED_TEST_SUITE(TileTest, DimTiles);

TYPED_TEST(TileTest, expectedNbrOfTilesPerDimToCoverTheBox)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(50)};
    auto const tile_size = PHARE::core::ConstArray<int, dim>(4);
    TypeParam cs{box, tile_size};


    auto const& shape = cs.shape();
    for (auto i = 0u; i < dim; ++i)
        EXPECT_EQ(shape[i], 13);
}

TYPED_TEST(TileTest, cluserSetSizeIsCorrect)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(47)};
    auto const tile_size = PHARE::core::ConstArray<int, dim>(4);
    TypeParam cs{box, tile_size};

    EXPECT_EQ(cs.size(), std::pow(12, dim));
}


TYPED_TEST(TileTest, totalTileSetSurfaceIsEqualToBoxSurface)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(49)};
    auto const tile_size = PHARE::core::ConstArray<int, dim>(4);
    TypeParam cs{box, tile_size};

    auto surface = 0.;
    for (auto i = 0u; i < cs.size(); ++i)
    {
        auto current_surface = 1.;
        for (auto d = 0u; d < dim; ++d)
        {
            auto l = (cs[i].upper[d] - cs[i].lower[d] + 1);
            current_surface *= l;
        }
        surface += current_surface;
    }
    EXPECT_EQ(surface, box.size());
}



TYPED_TEST(TileTest, tileHasNoOverlapWithOthers)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(54)};
    auto const tile_size = PHARE::core::ConstArray<int, dim>(4);
    TypeParam cs{box, tile_size};

    for (auto const& tile : cs)
    {
        for (auto const& other : cs)
        {
            if (&tile != &other)
            {
                auto const box1 = Box<int, dim>{tile.lower, tile.upper};
                auto const box2 = Box<int, dim>{other.lower, other.upper};
                auto overlap    = box1 * box2;
                EXPECT_FALSE(overlap.has_value());
            }
        }
    }
}


TYPED_TEST(TileTest, retrieveTilesFromBoxOverlap)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(54)};
    auto const tile_size = PHARE::core::ConstArray<int, dim>(4);
    TypeParam cs{box, tile_size};
    Box<int, dim> selection_box{ConstArray<int, dim>(11), ConstArray<int, dim>(34)};

    auto expected_nbr = std::pow(7, dim);
    auto overlapeds   = cs.export_overlaped_with(selection_box);
    EXPECT_EQ(overlapeds.size(), expected_nbr);

    auto completes   = 0.;
    auto incompletes = 0.;
    for (auto const& overlaped : overlapeds)
    {
        auto const& [is_complete, tile] = overlaped;
        if (is_complete)
            ++completes;
        else
            ++incompletes;
    }
    EXPECT_EQ(completes, std::pow(5, dim));
    EXPECT_EQ(incompletes, std::pow(7, dim) - std::pow(5, dim));
}


TYPED_TEST(TileTest, cannotCreateTileWithTileSizeBiggerThanBox)
{
    constexpr auto dim = TypeParam::dimension;
    Box<int, dim> box{ConstArray<int, dim>(0), ConstArray<int, dim>(5)};
    auto const tile_size = PHARE::core::ConstArray<int, dim>(7); // larger than box shape
    EXPECT_THROW(std::make_unique<TypeParam>(box, tile_size), std::runtime_error);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
