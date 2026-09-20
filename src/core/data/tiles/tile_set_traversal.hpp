#ifndef PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_TRAVERSAL_HPP

#include "core/utilities/box/box.hpp"

#include <array>
#include <cassert>
#include <optional>
#include <functional>
#include <type_traits>

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


// per-dimension "backward/mixed-sign" neighbour offsets: the (3**dim - 1) neighbours not
// already reachable via the (2**dim - 1) forward `_links` built by TileSet::build_links.
template<std::size_t dim>
auto constexpr tile_neighbour_backlinks();

template<>
auto constexpr tile_neighbour_backlinks<1>()
{
    return std::array<std::array<int, 1>, 1>{{
        {-1},
    }};
}

template<>
auto constexpr tile_neighbour_backlinks<2>()
{
    return std::array<std::array<int, 2>, 5>{{
        {-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {1, -1},
    }};
}

template<>
auto constexpr tile_neighbour_backlinks<3>()
{
    // 19 = (3**3) - 1 - 7
    return std::array<std::array<int, 3>, 19>{{
        {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0},
        {-1, 0, 1},   {-1, 1, -1}, {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1},
        {0, -1, 0},   {0, -1, 1},  {0, 0, -1},  {0, 1, -1},  {1, -1, -1},
        {1, -1, 0},   {1, -1, 1},  {1, 0, -1},  {1, 1, -1},
    }};
}

// visits every tile whose relative position to `tile` is one of the (3**dim - 1) neighbours,
// combining the forward `_links` (built by TileSet::build_links) with the remaining
// (mixed-sign) directions, looked up directly via the tileset's cell map.
template<typename TileSet_t, typename Tile_t, typename Fn>
void traverse_tile_neighbours(TileSet_t& tileset, Tile_t& tile, Fn fn)
{
    auto constexpr static dim = TileSet_t::dimension;

    auto const neighbour_point = [&](auto point) -> std::optional<Point<std::uint32_t, dim>> {
        for (std::uint8_t i = 0; i < dim; ++i)
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

    for (auto const& link : tile_neighbour_backlinks<dim>())
        if (auto const point = neighbour_point(Point{link}))
            fn(*tileset.at(*point));

    for (auto* link : tile.links())
        if (link)
            fn(*link);
}

// visits every ordered neighbouring tile pair (t0, t1) in the tileset exactly once per
// direction, i.e. both (A, B) and (B, A) are visited for adjacent tiles A and B.
template<typename TileSet_t>
void visit_tile_neighbours(TileSet_t& ts, auto fn)
{
    // link()/links() are declared on the (possibly base) linked type - e.g. GridTile owns its
    // tiles but only ever exposes FieldTile* through link(), so bind to that, not
    // TileSet_t::value_type, or a GridTileSet's neighbours won't convert back to GridTile&.
    using Tile   = std::remove_pointer_t<std::remove_reference_t<decltype(ts[0].link(0))>>;
    using TileFn = std::function<void(Tile&)>;

    TileFn const doX = [&](auto& t0) {
        traverse_tile_neighbours(ts, t0, [&](auto& t1) { fn(t0, t1); });

        if (auto l0 = t0.link(0))
            doX(*l0);
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
