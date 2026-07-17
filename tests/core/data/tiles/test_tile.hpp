
#ifndef PHARE_CORE_DATA_TILES_TEST_TILE_HPP
#define PHARE_CORE_DATA_TILES_TEST_TILE_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"


#include "gtest/gtest.h"

#include <unordered_set>


using namespace PHARE::core;


template<typename TileSet>
class TileTestBase
{
public:
    using TileSet_t = TileSet;

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
class TileTest : public TileTestBase<TileSet>, public ::testing::Test
{
public:
    TileTest()
        : TileTestBase<TileSet>{
              Box<int, TileSet::dimension>{ConstArray<int, TileSet::dimension>(0),
                                           ConstArray<int, TileSet::dimension>(47)},
              ConstArray<std::size_t, TileSet::dimension>(4)}
    {
    }
};


template<std::size_t dim>
class TileMock : public Box<int, dim>
{
    using Super = Box<std::int32_t, dim>;

public:
    TileMock(Super const& box)
        : Super{box}
    {
    }

    auto& link(std::size_t const idx) { return _links[idx]; }
    auto& links() { return _links; }
    auto& links() const { return _links; }

    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }

    std::array<TileMock*, 7> _links{ConstArray<TileMock*, 7>(nullptr)}; // 7 in 3d
};

#endif /*PHARE_CORE_DATA_TILES_TEST_TILE_HPP*/
