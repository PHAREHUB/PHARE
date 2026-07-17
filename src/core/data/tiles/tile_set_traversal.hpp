#ifndef PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP

#include "core/data/tiles/tile_set.hpp"
#include "core/utilities/box/box.hpp"

#include <cassert>
#include <cstdint>
#include <optional>

namespace PHARE::core
{

template<typename TileSet_t, typename Box_t, typename Fn>
void traverse_tiles(TileSet_t& tileset, Box_t const& box, Fn fn)
{
    using Tile_t = typename TileSet_t::value_type;

    std::function<void(Tile_t&)> const doX = [&](auto& tile) {
        fn(tile);

        if (auto nextX = tile.link(0); nextX and (box * (**nextX)))
            doX(*nextX);
    };
    std::function<void(Tile_t&)> const doY = [&](auto& tile) {
        doX(tile);

        if constexpr (TileSet_t::dimension > 1)
            if (auto nextY = tile.link(1); nextY and (box * (**nextY)))
                doY(*nextY);
    };
    std::function<void(Tile_t&)> const doZ = [&](auto& tile) {
        doY(tile);

        if constexpr (TileSet_t::dimension == 3)
            if (auto nextZ = tile.link(3); nextZ and (box * (**nextZ)))
                doZ(*nextZ);
    };

    doZ(*tileset.at(box.lower));
}


template<typename TileSet0, typename TileSet1, typename Box_t, typename Fn>
void traverse_tilesets_overlap(TileSet0& ts0, TileSet1& ts1, Box_t const& box, Fn fn)
{
    // box is in ts0 amr space
    using Tile0  = decltype(*ts0.views().at(box.lower));
    using Tile1  = decltype(*ts1().at(box.lower));
    using TileFn = std::function<void(Tile0&, Tile1&)>;

    TileFn const doX = [&](auto& t0, auto& t1) {
        if (box * *t0 and box * *t1 and *t0 * *t1)
            fn(t0, t1);

        auto const nt0 = t0.link(0);
        auto const nt1 = t1.link(0);

        if (nt0) // do y and z same
            doX(*nt0, t1);

        if (nt0 and nt1)
            doX(*nt0, *nt1);
    };
    TileFn const doY = [&](auto& t0, auto& t1) {
        doX(t0, t1);

        if constexpr (TileSet0::dimension > 1)
        {
            auto const nt0 = t0.link(1);
            auto const nt1 = t1.link(1);

            if (nt0 and nt1 and box * *nt0 and box * *nt1)
                doY(*nt0, *nt1);
        }
    };
    TileFn const doZ = [&](auto& t0, auto& t1) {
        doY(t0, t1);

        if constexpr (TileSet0::dimension == 3)
        {
            auto const nt0 = t0.link(3);
            auto const nt1 = t1.link(3);

            if (nt0 and nt1 and box * *nt0 and box * *nt1)
                doZ(*nt0, *nt1);
        }
    };

    doZ(*ts0.views().at(box.lower), *ts1().at(box.lower));
}


template<typename TileSet0, typename TileSet1, typename Box_t, typename Fn, typename Shift>
void traverse_tilesets_overlap(TileSet0& ts0, TileSet1& ts1, Box_t const& box, Fn fn, Shift shifter)
{
    // box is in ts0 amr space
    using Tile0  = decltype(*ts0.views().at(box.lower));
    using Tile1  = decltype(*ts1().at(box.lower));
    using TileFn = std::function<void(Tile0&, Tile1&)>;

    auto const ts1_box = shift(box, shifter * -1); // untransformed box in ts1 space

    TileFn const doX = [&](auto& t0, auto& t1) {
        if (box * *t0 and ts1_box * *t1)
            fn(t0, t1);

        auto const nt0 = t0.link(0);
        auto const nt1 = t1.link(0);

        if (nt0) // do y and z same
            doX(*nt0, t1);

        if (nt0 and nt1 and box * *nt0 and ts1_box * *nt1)
            doX(*nt0, *nt1);
    };
    TileFn const doY = [&](auto& t0, auto& t1) {
        doX(t0, t1);

        if constexpr (TileSet0::dimension > 1)
        {
            auto const nt0 = t0.link(1);
            auto const nt1 = t1.link(1);
            if (nt0 and nt1 and box * *nt0 and ts1_box * *nt1)
                doY(*nt0, *nt1);
        }
    };
    TileFn const doZ = [&](auto& t0, auto& t1) {
        doY(t0, t1);

        if constexpr (TileSet0::dimension == 3)
        {
            auto const nt0 = t0.link(3);
            auto const nt1 = t1.link(3);
            if (nt0 and nt1 and box * *nt0 and ts1_box * *nt1)
                doZ(*nt0, *nt1);
        }
    };

    doZ(*ts0.views().at(box.lower), *ts1().at(box.lower - shifter));
}


template<typename TileSet_t, typename Fn>
void traverse_ghost_boundary_tiles(TileSet_t& tileset, Fn fn)
{
    auto const domain_tile_box = shrink(tileset.box(), 2); // bleh

    for (auto& tile : tileset)
        if (auto overlap = domain_tile_box * tile; overlap != tile)
            fn(tile);
}

template<typename TileSet_t, typename Tile_t, typename Fn>
void traverse_tile_neighbours(TileSet_t& tileset, Tile_t& tile, Fn fn)
{
    if constexpr (TileSet_t::dimension == 3)
    {
        // 19 = (3**3) - 1 - 7
        auto constexpr backlinks = std::array<std::array<int, 3>, 19>{{
            {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0},
            {-1, 0, 1},   {-1, 1, -1}, {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1},
            {0, -1, 0},   {0, -1, 1},  {0, 0, -1},  {0, 1, -1},  {1, -1, -1},
            {1, -1, 0},   {1, -1, 1},  {1, 0, -1},  {1, 1, -1},
        }};

        auto const neighbour_point = [&](auto point) -> std::optional<Point<std::uint32_t, 3>> {
            for (std::uint8_t i = 0; i < 3; ++i)
            {
                if (point[i] < 0)
                    point[i] = tile.lower[i] - 1;
                else if (point[i] > 0)
                    point[i] = tile.upper[i] + 1;
                else
                    point[i] = tile.lower[i];
            }

            if (!isIn(point, tileset.box()))
                return std::nullopt;

            return (point - tileset.box().lower).as_unsigned();
        };

        for (auto const& link : backlinks)
            if (auto const point = neighbour_point(Point{link}))
                fn(*tileset.at(*point));

        for (auto* link : tile.links())
            if (link)
                fn(*link);
    }
    else
        throw std::runtime_error("TODO");
}

template<typename TileSet_t>
void visit_tile_neighbours(TileSet_t& ts, auto fn)
{
    using Tile   = TileSet_t::value_type;
    using TileFn = std::function<void(Tile&)>;

    if constexpr (TileSet_t::dimension < 3)
        throw std::runtime_error("TODO");

    TileFn const doX = [&](auto& t0) {
        traverse_tile_neighbours(ts, t0, [&](auto& t1) { fn(t0, t1); });

        // for (std::uint8_t i = 1; i < t0.links().size(); ++i)
        //     if (auto li = t0.link(i))
        //         fn(t0, *li);
        if (auto l0 = t0.link(0))
        {
            // fn(t0, *l0);
            doX(*l0);
        }
    };
    TileFn const doY = [&](auto& t0) {
        doX(t0);

        if constexpr (TileSet_t::dimension > 1)
            if (auto l1 = t0.link(1))
                doY(*l1);
    };
    TileFn const doZ = [&](auto& t0) {
        doY(t0);

        if constexpr (TileSet_t::dimension == 3)
            if (auto l3 = t0.link(3))
                doZ(*l3);
    };


    doZ(ts[0]);
}



} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP*/
