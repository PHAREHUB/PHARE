#include <ctype.h>
#include <string>

#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"
#include "core/data/tiles/tile_set.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>


using namespace PHARE::core;


template<typename TileSet>
class TileTestBase
{
public:
    static auto constexpr dimension = TileSet::dimension;

    TileTestBase(Box<int, dimension> box_, std::array<std::size_t, dimension> const& tile_size)
        : box{box_}
        , tileSet{box, tile_size}
    {
    }

    Box<int, dimension> box;
    TileSet tileSet;
};

template<typename TileSet>
class TileTestBoxShapeNotMultipleTileSize : public TileTestBase<TileSet>, public ::testing::Test
{
public:
    TileTestBoxShapeNotMultipleTileSize()
        : TileTestBase<TileSet>{
              Box<int, TileSet::dimension>{ConstArray<int, TileSet::dimension>(0),
                                           ConstArray<int, TileSet::dimension>(54)},
              ConstArray<std::size_t, TileSet::dimension>(4)}
    {
    }
};

template<typename TileSet>
class TileTestBoxShapeMultipleTileSize : public TileTestBase<TileSet>, public ::testing::Test
{
public:
    TileTestBoxShapeMultipleTileSize()
        : TileTestBase<TileSet>{
              Box<int, TileSet::dimension>{ConstArray<int, TileSet::dimension>(0),
                                           ConstArray<int, TileSet::dimension>(47)},
              ConstArray<std::size_t, TileSet::dimension>(4)}
    {
    }
};


template<typename TileSet>
class TileTest : public ::testing::Test
{
};


template<std::size_t dim>
class TileMock : public Box<int, dim>
{
};

using DimTiles = testing::Types<TileSet<TileMock<1>>, TileSet<TileMock<2>>, TileSet<TileMock<3>>>;

TYPED_TEST_SUITE(TileTestBoxShapeNotMultipleTileSize, DimTiles);
TYPED_TEST_SUITE(TileTestBoxShapeMultipleTileSize, DimTiles);
TYPED_TEST_SUITE(TileTest, DimTiles);


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, expectedNbrOfTilesPerDimToCoverTheBox)
{
    auto const& shape = this->tileSet.shape();
    for (auto i = 0u; i < this->dimension; ++i)
        EXPECT_EQ(shape[i], 14);
}

TYPED_TEST(TileTestBoxShapeMultipleTileSize, cluserSetSizeIsCorrect)
{
    EXPECT_EQ(this->tileSet.size(), std::pow(12, this->dimension));
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, totalTileSetSurfaceIsEqualToBoxSurface)
{
    auto surface = 0.;
    for (auto i = 0u; i < this->tileSet.size(); ++i)
    {
        auto current_surface = 1.;
        for (auto d = 0u; d < this->dimension; ++d)
        {
            auto l = (this->tileSet[i].upper[d] - this->tileSet[i].lower[d] + 1);
            current_surface *= l;
        }
        surface += current_surface;
    }
    EXPECT_EQ(surface, this->box.size());
}



TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, tileHasNoOverlapWithOthers)
{
    auto constexpr dim = TypeParam::dimension;
    for (auto const& tile : this->tileSet)
    {
        for (auto const& other : this->tileSet)
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


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, retrieveTilesFromBoxOverlap)
{
    auto constexpr dim = TypeParam::dimension;
    Box<int, dim> selection_box{ConstArray<int, dim>(11), ConstArray<int, dim>(34)};

    auto expected_nbr = std::pow(7, this->dimension);
    auto overlapeds   = this->tileSet.export_overlaped_with(selection_box);
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
    auto const tile_size = PHARE::core::ConstArray<std::size_t, dim>(7); // larger than box shape
    EXPECT_THROW(std::make_unique<TypeParam>(box, tile_size), std::runtime_error);
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, canRetrieveTileFromCell)
{
    auto constexpr dim = TypeParam::dimension;
    auto tile          = [&]() {
        if constexpr (dim == 1)
            return this->tileSet.at(13);
        else if constexpr (dim == 2)
            return this->tileSet.at(13, 13);
        else if constexpr (dim == 3)
            return this->tileSet.at(13, 13, 13);
    }();
    auto const expected_box = Box<int, dim>{ConstArray<int, dim>(12), ConstArray<int, dim>(15)};
    EXPECT_TRUE(*tile == expected_box);
}


TYPED_TEST(TileTestBoxShapeNotMultipleTileSize, getView)
{
    auto view = this->tileSet.make_view();
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
