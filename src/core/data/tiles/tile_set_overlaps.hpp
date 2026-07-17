#ifndef PHARE_CORE_DATA_TILES_TILE_SET_OVERLAPS_HPP
#define PHARE_CORE_DATA_TILES_TILE_SET_OVERLAPS_HPP


#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include <type_traits>
#include <unordered_set>


namespace PHARE::core
{
template<std::size_t dim, typename T>
struct NdSpanSet
{
    NdSpanSet(Box<std::uint32_t, dim> const box, std::size_t const size)
        : box_{box}
        , vec(size, 0)
        , displs(size, 0)

    {
    }

    auto& box() const { return box_; }

    Box<std::uint32_t, dim> box_;
    std::vector<T> vec;
    std::vector<default_span_size_t> displs;
    NdArrayVector<dim, Span<T>> cells{box_.shape()};
};

template<std::size_t dim>
struct TileBoxSpanSet
{
    using Tile_t = Box<int, dim>;

    TileBoxSpanSet(auto const& layout, auto gb, std::size_t const size = 0)
        : box_{layout.AMRBox()}
        , ghost_box_{gb}
        , vec(size, 0)
        , displs(size, 0)
    {
    }

    auto& box() const { return ghost_box_; }
    void resize(std::size_t const size)
    {
        vec.resize(size, 0);
        displs.resize(size, 0);
    }

    Box<int, dim> box_;
    Box<std::uint32_t, dim> ghost_box_;
    TileSet<Box<int, dim>> tiles{box_};
    std::vector<Tile_t*> vec;
    std::vector<default_span_size_t> displs;
    NdArrayVector<dim, Span<Tile_t*>> cells{ghost_box_.shape()};
};



template<typename GridLayout_t>
auto make_nd_span_set_from(GridLayout_t const& layout, auto const ghostboxer)
{
    auto constexpr static dim = GridLayout_t::dimension;

    using NdArrViewer = NdArrayViewer<dim>;

    TileBoxSpanSet<dim> spanset{layout, ghostboxer(layout)};

    auto const patch_ghost_box = ghostboxer(layout);

    NdArrayVector<dim, std::size_t> tiles_per_cell{patch_ghost_box.shape()};

    auto const map = [&](auto const fn) {
        for (auto& tile : spanset.tiles)
            for (auto const& bix : ghostboxer(layout.copy_as(tile)))
                fn(tile, bix);
    };

    map([&](auto&, auto const& bix) { tiles_per_cell(bix) += 1; });

    spanset.resize(sum(tiles_per_cell));

    std::size_t offset = 0;
    for (auto const& bix : patch_ghost_box)
    {
        spanset.cells(bix).s = tiles_per_cell(bix);
        auto const idx       = NdArrViewer::idx(patch_ghost_box.shape(), *bix);
        spanset.displs[idx]  = offset;
        offset += spanset.cells(bix).s;
    }

    tiles_per_cell.zero();

    map([&](auto& tile, auto const& bix) {
        auto const idx = NdArrViewer::idx(patch_ghost_box.shape(), *bix);
        spanset.vec[spanset.displs[idx] + tiles_per_cell(bix)] = &tile;
        ++tiles_per_cell(bix);
    });

    for (auto const& bix : patch_ghost_box)
        spanset.cells(bix).ptr
            = &spanset.vec[spanset.displs[NdArrViewer::idx(patch_ghost_box.shape(), *bix)]];

    return spanset;
}



template<typename GridLayout_t>
auto make_nd_span_set_for_qty(GridLayout_t const& layout, auto const pq)
{
    return make_nd_span_set_from(layout,
                                 [&](auto const& layout) { return layout.ghostBoxFor(pq); });
}


template<typename TileSet_t, typename GridLayout_t>
auto make_nd_span_set_from(TileSet_t& tiles, GridLayout_t const& layout, auto&& ghostboxer)
{
    static_assert(not std::is_const_v<std::remove_reference_t<decltype(tiles)>>);
    auto constexpr static dim      = GridLayout_t::dimension;
    auto constexpr static is_const = std::is_const_v<decltype(tiles)>;
    using value_type               = TileSet_t::value_type;
    using Tile_t      = value_type; // std::conditional_t<is_const, value_type const, value_type>;
    using NdArrViewer = NdArrayViewer<dim>;

    static_assert(not std::is_const_v<Tile_t>);
    assert(tiles.size());

    auto const amr_ghost_box = ghostboxer(layout.AMRBox());

    auto const lcl_box
        = [&](auto const& box) { return as_unsigned(shift(box, amr_ghost_box.lower * -1)); };

    auto const patch_ghost_box = lcl_box(amr_ghost_box);

    NdArrayVector<dim, std::size_t> tiles_per_cell{patch_ghost_box.shape()};

    auto map = [&](auto&& fn) {
        static_assert(not std::is_const_v<decltype(tiles)>);
        for (Tile_t& tile : tiles)
            for (auto const& bix : lcl_box(tile))
                fn(&tile, bix);
        for (Tile_t& tile : tiles)
            for (auto& box : ghostboxer(tile).remove(tile))
                for (auto const& bix : lcl_box(box))
                    fn(&tile, bix);
    };

    map([&](auto, auto const& bix) { tiles_per_cell(bix) += 1; });

    NdSpanSet<dim, Tile_t*> spanset{patch_ghost_box, sum(tiles_per_cell)};

    std::size_t offset = 0;
    for (auto const& bix : patch_ghost_box)
    {
        spanset.cells(bix).s = tiles_per_cell(bix);
        auto const idx       = NdArrViewer::idx(patch_ghost_box.shape(), *bix);
        spanset.displs[idx]  = offset;
        offset += spanset.cells(bix).s;
    }

    tiles_per_cell.zero();

    map([&](Tile_t* tile, auto const& bix) {
        auto const idx = NdArrViewer::idx(patch_ghost_box.shape(), *bix);
        spanset.vec[spanset.displs[idx] + tiles_per_cell(bix)] = tile;
        ++tiles_per_cell(bix);
    });

    for (auto const& bix : patch_ghost_box)
        spanset.cells(bix).ptr
            = &spanset.vec[spanset.displs[NdArrViewer::idx(patch_ghost_box.shape(), *bix)]];

    return spanset;
}

template<auto quantity, typename TileSet_t, typename GridLayout_t>
auto make_qty_nd_span_set_from(TileSet_t& tiles, GridLayout_t const& layout)
{
    return make_nd_span_set_from(tiles, layout, [&](auto const& box) {
        return layout.copy_as(box).AMRGhostBoxFor(quantity);
    });
}

template<typename TileSet_t, typename GridLayout_t>
auto make_qty_nd_span_set_from(TileSet_t& tiles, GridLayout_t const& layout, auto const qty)
{
    return make_nd_span_set_from(
        tiles, layout, [&](auto const& box) { return layout.copy_as(box).AMRGhostBoxFor(qty); });
}


template<typename TileSet_t, typename GridLayout_t>
auto make_particle_nd_span_set_from(TileSet_t& tiles, GridLayout_t const& layout)
{
    return make_nd_span_set_from(tiles, layout, [](auto const& box) {
        return grow(box, GridLayout_t::nbrParticleGhosts());
    });
}


template<std::size_t dim, typename T>
auto unique_tiles_for(NdSpanSet<dim, T> const& spanset, Box<std::uint32_t, dim> const& box)
{
    assert(spanset.box() * box);
    std::unordered_set<T> tiles;
    for (auto const& bix : box)
        for (auto& tile_ptr : spanset.cells(bix))
            tiles.insert(tile_ptr);

    return tiles;
}


template<std::size_t dim, typename T>
auto unique_tiles_for(NdSpanSet<dim, T> const& spanset, Point<std::uint32_t, dim> const& point)
{
    return unique_tiles_for(spanset, asBox(point));
}



template<std::size_t dim, typename T>
void onBox(NdSpanSet<dim, T> const& spanset, Box<std::uint32_t, dim> const& box, auto&& fn)
{
    for (auto& tp : unique_tiles_for(spanset, box))
        fn(*tp);
}


void check_tile_span_set(auto const& spanset)
{
    for (auto const& span : spanset.cells)
        if (span.size() == 0)
            throw std::runtime_error("bad span");

    for (auto const& span : spanset.cells)
        for (auto const& p : span)
            if (!p)
                throw std::runtime_error("bad pointer");
}


} // namespace PHARE::core


#endif
