#ifndef PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_MAPPER_HPP

#include "core/utilities/point/point.hpp"
#include "core/data/tiles/tiling_options.hpp"

#include <vector>
#include <cassert>
#include <algorithm>

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
    template<typename... Args>
    auto map(Args&&... args);

    TilingOptions const opts = TilingOptions::from_env();
};


template<typename TileSet_t>
template<typename... Args>
auto Tiler<TileSet_t>::map(Args&&... args)
{
    auto constexpr static dim = TileSet_t::dimension;
    using Box_t               = TileSet_t::Box_t;
    auto const& shape         = this->box.shape();
    auto& tiles               = this->tile_set();

    // Greedily fills with max size chunks, sizes are multiples of `step` (1 for level 0, 2 for
    // fine levels). A remainder below min size is rebalanced with the previous tile, so no tile
    // ever exceeds max size or is empty, even if min/max don't round to step (e.g. min=max=15)
    auto const split_1d = [&](std::size_t const length,
                              std::size_t const step) -> std::vector<std::size_t> {
        assert(length % step == 0);
        if (length < opts.min_patch_size_before_split)
            return {length};


        auto const round_up        = [&](auto const v) { return (v + step - 1) / step * step; };
        std::size_t const max_size = std::max(step, opts.max_tile_size / step * step);
        std::size_t const min_size = std::min(max_size, round_up(opts.min_tile_size));
        std::vector<std::size_t> sizes;
        std::size_t rem = length;
        while (rem > max_size)
        {
            sizes.push_back(max_size);
            rem -= max_size;
        }
        // rem is a multiple of step, in [step, max_size]
        if (rem >= min_size or sizes.empty())
            sizes.push_back(rem);
        else
        {
            auto const total = sizes.back() + rem; // in (max_size, 2 * max_size)
            sizes.back()     = round_up(total / 2);
            sizes.push_back(total - sizes.back());
        }
        return sizes;
    };

    // Even box shape → fine level → even tiles required, so that each tile's outermost ghost
    // cell lands on an even AMR index and is therefore filled by the init refiner;
    // odd shape → level 0 → any tiles.
    bool const all_even_shape = [&] {
        for (std::size_t d = 0; d < dim; ++d)
            if (shape[d] % 2 != 0)
                return false;
        return true;
    }();

    auto const do_split = [&](std::size_t const length) {
        return split_1d(length, all_even_shape ? 2 : 1);
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
void tile_set_make_tiles(TileSet_t& tile_set, TilingOptions const& opts, Args&&... args)
{
    Tiler<TileSet_t>{{tile_set}, opts.validate()}.map(args...);
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
