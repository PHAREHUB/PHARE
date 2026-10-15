

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set_mapper.hpp"

#include "test_tile.hpp"

#include <gtest/gtest.h>


template<typename TileSet_>
struct TileSetInfo
{
    auto static constexpr dimension = TileSet_::dimension;
    using Tile_t                    = TileSet_::value_type;
    using TileSet_t                 = TileSet_;

    Box<int, TileSet_t::dimension> box;
    std::array<std::uint32_t, TileSet_t::dimension> shape;
    std::array<std::size_t, TileSet_t::dimension> tile_size;
};

template<typename TileSet_>
class TileMappingTest : public ::testing::Test
{
public:
    auto static constexpr dimension = TileSet_::dimension;
    using Tile_t                    = TileSet_::value_type;
    using TileSet_t                 = TileSet_;
    using info_t                    = TileSetInfo<TileSet_t>;

    TileMappingTest() {}

    auto static count_cells(TileSet_t const& tiles)
    {
        return sum_from(tiles, [](auto const& e) { return e.size(); });
    }

    auto static all_in(TileSet_t const& tiles, Box<int, dimension> const& box)
    {
        for (auto const& tile : tiles)
            if (box * tile != tile)
                return false;
        return true;
    };
};

template<typename Tile>
struct TestTileSet
{
    auto static constexpr dimension = Tile::dimension;
    using value_type                = Tile;
    using Tile_t                    = Tile;
    using Box_t                     = Box<int, dimension>;

    NO_DISCARD auto& box() const { return box_; }
    NO_DISCARD auto shape() const { return shape_; }
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return tiles_.size(); }
    NO_DISCARD auto begin() { return tiles_.begin(); }
    NO_DISCARD auto begin() const { return tiles_.begin(); }
    NO_DISCARD auto end() { return tiles_.end(); }
    NO_DISCARD auto end() const { return tiles_.end(); }
    NO_DISCARD auto& operator()() { return tiles_; }
    NO_DISCARD auto& operator()() const { return tiles_; }
    NO_DISCARD auto& operator[](std::size_t const i) _PHARE_ALL_FN_ { return tiles_[i]; }
    NO_DISCARD auto& operator[](std::size_t const i) const _PHARE_ALL_FN_ { return tiles_[i]; }

    Box<int, dimension> box_;
    std::array<std::uint32_t, dimension> shape_;
    std::array<std::size_t, dimension> tile_size;
    std::vector<Tile> tiles_;
};

using DimTiles = testing::Types<TestTileSet<Box<int, 3>>>;

TYPED_TEST_SUITE(TileMappingTest, DimTiles);

TYPED_TEST(TileMappingTest, view_outputs)
{
    auto static constexpr dimension = TestFixture::dimension;
    using TileSet_t                 = TestFixture::TileSet_t;
    using Box_t                     = TileSet_t::Box_t;

    auto const doBox = [&](auto const box, TilingOptions const& opts) {
        TileSet_t tiles{box};
        Tiler<TileSet_t>{{tiles}, opts}.map();

        EXPECT_GT(tiles.size(), 0);

        bool const even = for_N_all<dimension>([&](auto i) { return box.shape()[i] % 2 == 0; });
        for (auto const& tile : tiles)
            for (std::size_t i = 0; i < dimension; ++i)
            {
                auto const s = static_cast<std::size_t>(tile.shape()[i]);
                EXPECT_GT(s, 0u);
                if (static_cast<std::size_t>(box.shape()[i]) >= opts.min_patch_size_before_split)
                    EXPECT_LE(s, opts.max_tile_size);
                if (even)
                    EXPECT_EQ(s % 2, 0u);
            }

        EXPECT_TRUE(this->all_in(tiles, box));
        EXPECT_FALSE(any_overlaps_in(tiles, [](auto const& tile) { return tile; }));
        EXPECT_EQ(this->count_cells(tiles), product(box.shape()));
    };

    auto const options = std::vector<TilingOptions>{
        {4, 6, 8}, {4, 4, 8}, {4, 10, 8}, {5, 7, 10}, {15, 15, 30}, {6, 15, 12}};

    for (auto const& opts : options)
        for (std::uint8_t i = 4; i < 33; ++i)
            for (std::uint8_t j = 4; j < 33; ++j)
                for (std::uint8_t k = 4; k < 33; ++k)
                    doBox(Box_t{ConstArray<int, dimension>(0), {i - 1, j - 1, k - 1}}, opts);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
