#ifndef PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP

#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include <cassert>
#include <cstdint>

namespace PHARE::core
{

template<typename TileSet_t>
struct ATiler
{
    using Tile_t                    = TileSet_t::value_type;
    using Box_t                     = TileSet_t::Box_t;
    static auto constexpr dimension = TileSet_t::dimension;

    TileSet_t& tile_set;
    Box_t const& box = tile_set.box();
};


template<typename TileSet_t>
struct Tiler : public ATiler<TileSet_t>
{
    static inline std::string const min_key = "PHARE_TILING_MIN_BEFORE_SPLIT";
    static inline std::string const max_key = "PHARE_TILING_MAX_SIZE";

    template<typename... Args>
    auto f(Args&&... args);

    std::size_t const min_tile_size = 4;
    std::size_t const max_tile_size = get_env_as(max_key, 6); // matters for gpu dynamic shared mem
    std::size_t const min_patch_size_before_split = get_env_as(min_key, min_tile_size * 2);
};


template<typename TileSet_t>
template<typename... Args>
auto Tiler<TileSet_t>::f(Args&&... args)
{
    auto constexpr static dim          = TileSet_t::dimension;
    using Box_t                        = TileSet_t::Box_t;
    auto const& shape                  = this->box.shape();
    auto& tiles                        = this->tile_set();
    std::size_t const min_before_split = min_patch_size_before_split;

    // Returns tile sizes for an arbitrary-length dimension (may produce odd sizes).
    // Used for level 0 patches that can have any cell count.
    auto split_1d = [&](std::size_t const length) {
        if (length < min_before_split)
            return std::vector<std::size_t>{length};
        if (length == 8)
            return std::vector<std::size_t>{4, 4};
        if (length == 9)
            return std::vector<std::size_t>{5, 4};
        if (length == 10)
            return std::vector<std::size_t>{5, 5};
        if (length < 13)
            return std::vector<std::size_t>{6, length - 6};

        std::size_t const border = 4;
        std::size_t middle       = length - (2 * border);
        std::vector<std::size_t> sizes{border};

        if (middle / 4 < 3)
            sizes.emplace_back(middle);
        else
        {
            if (middle / 6 < 3)
            {
                sizes.emplace_back(6);
                sizes.emplace_back(middle - 6);
            }
            else
            {
                auto rem = middle % 6;
                while (middle > 5)
                {
                    std::size_t add = rem > 1 ? 2 : rem > 0 ? 1 : 0;
                    sizes.emplace_back(6 + add);
                    middle -= (6 + add);
                    rem -= add;
                }
                assert(middle == 0);
            }
        }

        sizes.push_back(border);
        return sizes;
    };

    // Returns tile sizes where every size is even. Required for fine-level patches
    // (even AMR alignment) so that each tile's outermost ghost cell lands on an even
    // AMR index and is therefore filled by the init refiner.
    auto split_1d_even = [&](std::size_t const length) -> std::vector<std::size_t> {
        assert(length % 2 == 0);
        if (length < min_before_split)
            return {length};
        std::size_t const max_even = max_tile_size - (max_tile_size % 2);
        assert(max_even >= 6);
        std::vector<std::size_t> sizes;
        std::size_t rem = length;
        while (rem > max_even)
        {
            sizes.push_back(max_even);
            rem -= max_even;
        }
        // rem is even, in [2, max_even]
        if (rem >= 4)
            sizes.push_back(rem);
        else // rem == 2: borrow 2 from the previous tile to form a trailing tile of 4
        {
            assert(!sizes.empty());
            sizes.back() -= 2;
            sizes.push_back(4);
        }
        return sizes;
    };

    // Even box shape → fine level → even tiles required; odd shape → level 0 → any tiles.
    bool const all_even_shape = [&] {
        for (std::size_t d = 0; d < dim; ++d)
            if (shape[d] % 2 != 0)
                return false;
        return true;
    }();

    auto const do_split = [&](std::size_t const length) {
        return all_even_shape ? split_1d_even(length) : split_1d(length);
    };

    auto const get_ranges = [](auto const& sizes) {
        std::vector<std::pair<std::size_t, std::size_t>> ranges;
        ranges.reserve(sizes.size());
        std::size_t current = 0;
        for (auto const s : sizes)
        {
            ranges.emplace_back(current, current + s);
            current += s;
        }
        return ranges;
    };

    auto const& box = this->box;

    auto const subdivide_1d = [&](auto const size_x) {
        auto const x_ranges = get_ranges(do_split(size_x));
        tiles.reserve(x_ranges.size());
        for (auto const& [xlo, xhi] : x_ranges)
            tiles.emplace_back(Box_t{Point{xlo} + box.lower, Point{xhi} + box.lower - 1}, args...);
    };

    auto const subdivide_2d = [&](auto const size_x, auto const size_y) {
        auto const x_ranges = get_ranges(do_split(size_x));
        auto const y_ranges = get_ranges(do_split(size_y));
        tiles.reserve(x_ranges.size() * y_ranges.size());
        for (auto const& [xlo, xhi] : x_ranges)
            for (auto const& [ylo, yhi] : y_ranges)
                tiles.emplace_back(
                    Box_t{Point{xlo, ylo} + box.lower, Point{xhi, yhi} + box.lower - 1}, args...);
    };

    auto const subdivide_3d = [&](auto const size_x, auto const size_y, auto const size_z) {
        auto const x_ranges = get_ranges(do_split(size_x));
        auto const y_ranges = get_ranges(do_split(size_y));
        auto const z_ranges = get_ranges(do_split(size_z));
        tiles.reserve(x_ranges.size() * y_ranges.size() * z_ranges.size());
        std::size_t count = 0;
        for (auto const& [xlo, xhi] : x_ranges)
            for (auto const& [ylo, yhi] : y_ranges)
                for (auto const& [zlo, zhi] : z_ranges)
                    tiles.emplace_back(Box_t{Point{xlo, ylo, zlo} + box.lower,
                                             Point{xhi, yhi, zhi} + box.lower - 1},
                                       args...),
                        ++count;

        assert(tiles.capacity() == count);
    };

    if constexpr (dim == 1)
        subdivide_1d(shape[0]);
    if constexpr (dim == 2)
        subdivide_2d(shape[0], shape[1]);
    if constexpr (dim == 3)
        subdivide_3d(shape[0], shape[1], shape[2]);

    assert(
        !any_overlaps_in(tiles, [](auto const& tile) { return static_cast<Box<int, dim>>(tile); }));

    assert(tiles.size() > 0);
}

template<typename TileSet_t, typename... Args>
void tile_set_make_tiles(TileSet_t& tile_set, Args&&... args)
{
    Tiler<TileSet_t>{{tile_set}}.f(args...);
}

template<typename TileSet_t, typename TileSet0>
TileSet_t tile_set_make_from_tiles(TileSet_t const& from, auto&&... args)
{
    TileSet_t out;
    for (auto const& in : from())
        out.tiles.emplace_back(in, args...);
    return out;
}

} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP*/
