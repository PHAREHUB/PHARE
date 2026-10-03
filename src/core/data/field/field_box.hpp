#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_HPP

#include "core/data/grid/grid.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"

#include <vector>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace PHARE::core
{
template<typename D>
struct FieldBorderSumOp : public PlusEquals<D>
{
};

// Picks whichever side is "more complete" (domain deposit vs. a partial/
// level-ghost-approximated contribution). For non-negative quantities (density)
// the more complete side is also the numerically larger one, so a plain max
// works. For signed quantities (e.g. bulk velocity) that's not true, so we
// compare magnitudes instead and keep the winning side's sign.
template<typename D>
struct FieldBorderMaxOp
{
    using value_type = D;

    void operator()(auto const& d0)
    {
        if (std::abs(d0) > std::abs(d))
            d = d0;
    }

    D& d;
};

template<typename Op>
bool constexpr is_border_op()
{
    using value_type = Op::value_type;
    return is_any_of<Op, FieldBorderSumOp<value_type>, FieldBorderMaxOp<value_type>>();
}

template<typename Op>
auto constexpr static is_border_op_v = is_border_op<Op>();

template<typename Op>
bool constexpr is_border_max_op()
{
    using value_type = Op::value_type;
    return std::is_same_v<Op, FieldBorderMaxOp<value_type>>;
}

template<typename Op>
auto constexpr static is_border_max_op_v = is_border_max_op<Op>();

template<typename Field_t>
class FieldBox
{
    using value_type = std::decay_t<typename Field_t::type>;
    using SetEqualOp = core::SetEqual<value_type>;

public:
    auto constexpr static dimension = Field_t::dimension;

    FieldBox(Field_t& field_, Box<int, dimension> const& amr_box_,
             Box<std::uint32_t, dimension> const& lcl_box_)
        : field{field_}
        , amr_box{amr_box_}
        , lcl_box{lcl_box_}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout)
        : FieldBox{field_, layout.AMRBox(), layout.ghostBoxFor(field_)}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout,
             Box<std::uint32_t, dimension> const& selection)
        : FieldBox{field_, layout.AMRBox(), selection}
    {
    }

    template<typename GridLayout_t>
    FieldBox(Field_t& field_, GridLayout_t const& layout, Box<int, dimension> const& selection)
        : FieldBox{field_, layout.AMRBox(), layout.AMRToLocal(selection)}
    {
    }

    template<typename Operator = SetEqualOp, typename Field_t0>
    void op(FieldBox<Field_t0> const& that);

    template<typename Operator = SetEqualOp>
    void op(value_type const val);

    template<typename Operator = SetEqualOp>
    void append_to(std::vector<value_type>& vec) const;

    auto& offset(auto const& offset)
    {
        offset_ = offset;
        return *this;
    }

    Field_t& field;
    Box<int, dimension> amr_box;
    Box<std::uint32_t, dimension> lcl_box;
    Point<int, dimension> offset_ = ConstArray<int, dimension>();
};

template<typename Field_t>
FieldBox(Field_t&, auto const&) -> FieldBox<Field_t>;
template<typename Field_t>
FieldBox(Field_t&, auto const&, auto const&) -> FieldBox<Field_t>;

template<typename Box_t>
struct BoxExpander // expands border tiles to cover the patch ghost box
{
    auto constexpr static dimension = Box_t::dimension;

    auto operator()(auto const& inbox) const
    {
        auto box           = inbox;
        auto const per_dim = [&](auto const di) {
            if (inbox.lower[di] == patch_box.lower[di])
                box.lower[di] = ghost_box.lower[di];
            auto const hi = inbox.upper[di];
            if (inbox.upper[di] == patch_box.upper[di])
                box.upper[di] = ghost_box.upper[di];
        };
        for (std::uint8_t i = 0; i < dimension; ++i)
            per_dim(i);
        return box;
    }

    Box_t ghost_box, patch_box;
};

template<typename Box_t>
BoxExpander(Box_t, Box_t) -> BoxExpander<Box_t>;

template<typename Operator, typename GridLayout_t, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<GridLayout_t, Args...>>& dst,
                       FieldBox<GridTileSet<GridLayout_t, Args...> const> const& src);

template<typename Tiles>
auto reduce_single(Tiles const& input)
    requires(is_field_tile_set_v<Tiles>);


template<typename Operator, typename GridLayout_t, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<GridLayout_t, Args...>>& dst,
                       FieldBox<GridTileSet<GridLayout_t, Args...> const> const& src)
    requires(is_border_op_v<Operator>)
{
    PHARE_LOG_SCOPE(3, "operate_on_fields_border_sum<GridTileSet,GridTileSet>");

    static_assert(dependent_false_v<Operator>);

    auto const pq = dst.field.physicalQuantity();
    assert(src.field.physicalQuantity() == pq);

    auto const skip = [](auto const& field, auto const& tile) {
        auto const tbox = grow(tile, 1);
        return *(tbox * field.layout().AMRBox()) == tbox;
    };
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box
        = shift(dst.field.layout().localToAMR(dst.lcl_box), src.offset_ * -1);

    for (auto& dst_tile : dst.field())
        if (!skip(dst.field, dst_tile))
            if (auto const dst_overlap
                = dst_selection_box * shift(get_box(dst_tile), src.offset_ * -1))
                for (auto const& src_tile : src.field())
                    if (!skip(src.field, src_tile))
                        if (auto const src_overlap = src_selection_box * get_box(src_tile))
                            if (auto const overlap = *src_overlap * *dst_overlap)
                            {
                                auto const lcl_src_box = src_tile.layout().AMRToLocal(*overlap);
                                auto const lcl_dst_box
                                    = dst_tile.layout().AMRToLocal(shift(*overlap, src.offset_));

                                assert(lcl_src_box.size() <= src_tile().size());
                                assert(lcl_dst_box.size() <= dst_tile().size());

                                auto src_it = lcl_src_box.begin();
                                auto dst_it = lcl_dst_box.begin();
                                for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end();
                                     ++src_it, ++dst_it)
                                {
                                    Operator{dst_tile()(*dst_it)}(src_tile()(*src_it));
                                }
                            }
}


template<typename Operator, typename GridLayout_t, typename... Args>
void operate_on_fields(FieldBox<GridTileSet<GridLayout_t, Args...>>& dst,
                       FieldBox<GridTileSet<GridLayout_t, Args...> const> const& src)
{
    PHARE_LOG_SCOPE(3, "operate_on_fields<GridTileSet,GridTileSet>");

    auto const pq = dst.field.physicalQuantity();
    assert(src.field.physicalQuantity() == pq);

    auto get_dst_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    auto get_src_box = [&](auto const& tile) { return tile.ghost_box(); };

    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box
        = shift(dst.field.layout().localToAMR(dst.lcl_box), src.offset_ * -1);

    assert(src_selection_box.shape() == dst_selection_box.shape());

    for (auto& dst_tile : dst.field())
        if (auto const dst_overlap
            = dst_selection_box * shift(get_dst_box(dst_tile), src.offset_ * -1))
            for (auto const& src_tile : src.field())
                if (auto const src_overlap = src_selection_box * get_src_box(src_tile))
                    if (auto const overlap = *src_overlap * *dst_overlap)
                    {
                        auto const lcl_src_box = src_tile.layout().AMRToLocal(*overlap);
                        auto const lcl_dst_box
                            = dst_tile.layout().AMRToLocal(shift(*overlap, src.offset_));

                        assert(lcl_src_box.size() <= src_tile().size());
                        assert(lcl_dst_box.size() <= dst_tile().size());

                        auto src_it = lcl_src_box.begin();
                        auto dst_it = lcl_dst_box.begin();
                        for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end();
                             ++src_it, ++dst_it)
                        {
                            Operator{dst_tile()(*dst_it)}(src_tile()(*src_it));
                        }
                    }
}


template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<GridTileSet<T0s...>>& dst, FieldBox<T1s...> const& src)
{
    PHARE_LOG_SCOPE(3, "operate_on_fields<GridTileSet,T1s...>");

    using Src = std::decay_t<decltype(src.field)>;
    static_assert(is_field_v<Src>);

    auto const amr_selection_box = dst.field.layout().localToAMR(dst.lcl_box);
    auto const pq                = dst.field.physicalQuantity();

    for (auto& dst_tile : dst.field())
    {
        auto const dst_tile_ghost_box = dst_tile.layout().AMRGhostBoxFor(pq);
        if (auto const overlap = amr_selection_box * dst_tile_ghost_box)
        {
            // src.lcl_box covers the full intersection in src-local coords.
            // Offset by how far this tile's overlap starts past amr_selection_box.lower.
            auto lcl_src_lower    = src.lcl_box.lower;
            auto const amr_offset = overlap->lower - amr_selection_box.lower;
            for (std::size_t i = 0; i < decltype(amr_offset)::dimension; ++i)
                lcl_src_lower[i] += static_cast<std::uint32_t>(amr_offset[i]);
            auto const overlap_shape = overlap->shape().template as<std::uint32_t>();
            auto const lcl_src_box   = Box<std::uint32_t, decltype(amr_offset)::dimension>{
                lcl_src_lower, lcl_src_lower + overlap_shape - 1};
            auto const lcl_dst_box = dst_tile.layout().AMRToLocal(*overlap);
            auto src_it            = lcl_src_box.begin();
            auto dst_it            = lcl_dst_box.begin();
            for (; dst_it != lcl_dst_box.end(); ++src_it, ++dst_it)
                Operator{dst_tile()(*dst_it)}(src.field(*src_it));
        }
    }
}

template<typename Field_t>
template<typename Operator, typename Field_t0>
void FieldBox<Field_t>::op(FieldBox<Field_t0> const& that)
{
    operate_on_fields<Operator>(*this, that);
}


// Max-combine: source tile ghost regions hold redundant *complete* copies of a node's
// value (not partial contributions), so overlapping tiles must be de-duplicated into one
// flat, unambiguous view before combining. Iterating src.field()'s tiles directly and
// applying Operator once per overlapping tile applies it multiple times to the same dst
// cell whenever several tiles' redundant ghost copies overlap the same destination region
// (last tile in container order wins, an order not guaranteed stable between serial and
// parallel runs, and not guaranteed stable across different tile decompositions).
template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<Grid<T0s...>>& dst, FieldBox<GridTileSet<T1s...> const> const& src)
    requires(is_border_max_op_v<Operator>)
{
    PHARE_LOG_SCOPE(3, "operate_on_fields_border_max<Grid, GridTileSet>");

    // de-duplicate only the border-touching tiles into one flat, unambiguous view. This
    // overload is only ever used to fill a border/ghost selection box, which non-border
    // (fully interior) tiles can never overlap, so populating them would be wasted work.
    using Tiles_t = GridTileSet<T1s...>;
    using Grid_t0 = typename Tiles_t::grid_type;

    auto const pq           = src.field.physicalQuantity();
    auto const patch_layout = src.field.layout();
    auto const patch_box    = patch_layout.AMRBox();
    auto const ghost_box    = patch_layout.AMRGhostBoxFor(pq);
    auto const expander     = BoxExpander{ghost_box, patch_box};

    Grid_t0 srcScratch{src.field.name(), pq, src.field.shape()};
    srcScratch.zero();
    for (auto const& tile : src.field())
    {
        auto const tile_gb = tile.layout().AMRGhostBoxFor(pq);
        if (*(tile_gb * patch_box) == tile_gb) // fully interior, cannot reach a border box
            continue;
        for (auto const& bix : expander(*tile))
        {
            auto const lix      = (bix - ghost_box.lower).as_unsigned();
            auto const tile_lix = (bix - tile_gb.lower).as_unsigned();
            srcScratch(lix)     = tile()(tile_lix);
        }
    }

    auto const dst_layout = src.field.layout().copy_as(dst.amr_box);
    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box = shift(dst_layout.localToAMR(dst.lcl_box), src.offset_ * -1);
    assert(src_selection_box.shape() == dst_selection_box.shape());

    if (auto const overlap = src_selection_box * dst_selection_box)
    {
        auto const lcl_src_box = src.field.layout().AMRToLocal(*overlap);
        auto const lcl_dst_box = dst_layout.AMRToLocal(shift(*overlap, src.offset_));

        auto src_it = lcl_src_box.begin();
        auto dst_it = lcl_dst_box.begin();
        for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end(); ++src_it, ++dst_it)
            Operator{dst.field(*dst_it)}(srcScratch(*src_it));
    }
}

// Sum-combine: each source tile's ghost region holds a genuinely *partial* contribution
// (e.g. raw particle deposit split across tile-ghost boundaries), so every overlapping
// tile must be visited and accumulated.
template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<Grid<T0s...>>& dst, FieldBox<GridTileSet<T1s...> const> const& src)
    requires(is_border_op_v<Operator> and !is_border_max_op_v<Operator>)
{
    // USED IN BORDER SUM SCHEDULES FOR TILES!
    PHARE_LOG_SCOPE(3, "operate_on_fields_border_sum<Grid, GridTileSet>");

    auto const pq         = dst.field.physicalQuantity();
    auto const dst_layout = src.field.layout().copy_as(dst.amr_box);
    auto const skip       = [](auto const& field, auto const& tile) {
        auto const tbox = grow(tile, 1);
        return *(tbox * field.layout().AMRBox()) == tbox;
    };
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    auto const src_selection_box = src.field.layout().localToAMR(src.lcl_box);
    auto const dst_selection_box = shift(dst_layout.localToAMR(dst.lcl_box), src.offset_ * -1);
    assert(src_selection_box.shape() == dst_selection_box.shape());

    for (auto const& src_tile : src.field())
    {
        if (skip(src.field, src_tile))
            continue;
        auto const src_tile_ghost_box = get_box(src_tile);
        if (auto const src_overlap = src_selection_box * src_tile_ghost_box)
            if (auto const overlap = *src_overlap * dst_selection_box)
            {
                auto const lcl_src_box = src_tile.layout().AMRToLocal(*overlap);
                auto const lcl_dst_box = dst_layout.AMRToLocal(shift(*overlap, src.offset_));

                assert(lcl_src_box.size() <= src_tile().size());

                auto src_it = lcl_src_box.begin();
                auto dst_it = lcl_dst_box.begin();
                for (; dst_it != lcl_dst_box.end() and src_it != lcl_src_box.end();
                     ++src_it, ++dst_it)
                    Operator{dst.field(*dst_it)}(src_tile()(*src_it));
            }
    }
}

template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<Grid<T0s...>>& dst, FieldBox<GridTileSet<T1s...> const> const& src)
{
    throw std::runtime_error(
        "FieldBox<Grid<T0s...>>& dst, FieldBox<GridTileSet<T1s...> const> const& src");
}

template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<Grid<T0s...>>& dst, FieldBox<FieldTileSet<T1s...> const> const& src)
{
    throw std::runtime_error("operate_on_fields(FieldBox<Grid<T0s...>>& dst, "
                             "FieldBox<FieldTileSet<T1s...> const> const& src)");
}


// final fallthroughs, only supports fields without tiles
template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<T0s...>& dst, FieldBox<T1s...> const& src)
{
    auto src_it = src.lcl_box.begin();
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++src_it, ++dst_it)
        Operator{dst.field(*dst_it)}(src.field(*src_it));
}
template<typename Operator, typename... T0s, typename... T1s>
void operate_on_fields(FieldBox<T0s...>&& dst, FieldBox<T1s...> const& src)
{
    operate_on_fields<Operator>(dst, src);
}

template<typename Operator, typename... Args>
void set_on_fields(FieldBox<FieldTileSet<Args...>>& dst, auto const val)
{
    if (dst.field().size() == 0)
        return;

    auto const pq = dst.field.physicalQuantity();
    auto const amr_selection_box
        = dst.field()[0].layout().copy_as(dst.amr_box).localToAMR(dst.lcl_box);

    for (auto& dst_tile : dst.field())
    {
        auto const dst_tile_ghost_box = dst_tile.layout().AMRGhostBoxFor(pq);
        if (auto const overlap = amr_selection_box * dst_tile_ghost_box)
        {
            auto const lcl_dst_box = dst_tile.layout().AMRToLocal(*overlap);
            for (auto dst_it = lcl_dst_box.begin(); dst_it != lcl_dst_box.end(); ++dst_it)
                Operator{dst_tile()(*dst_it)}(val);
        }
    }
}

template<typename Operator, typename... Args>
void set_on_fields(FieldBox<Args...>& dst, auto const val)
{
    auto dst_it = dst.lcl_box.begin();
    for (; dst_it != dst.lcl_box.end(); ++dst_it)
        Operator{dst.field(*dst_it)}(val);
}

template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::op(value_type const val)
{
    set_on_fields<Operator>(*this, val);
}

template<typename Operator, typename... Args>
void append_from_fields(FieldBox<FieldTileSet<Args...> const> const& src, auto& vec)
{
    if (src.field().size() == 0)
        return;

    auto const pq = src.field.physicalQuantity();
    auto const amr_selection_box
        = src.field()[0].layout().copy_as(src.amr_box).localToAMR(src.lcl_box);

    auto const base = vec.size();
    vec.resize(base + amr_selection_box.size(), 0);
    auto view = make_array_view(vec.data() + base, *amr_selection_box.shape().as_unsigned());

    using value_type = std::decay_t<decltype(view(amr_selection_box.lower.as_unsigned()))>;
    if constexpr (std::is_base_of_v<PlusEquals<value_type>, Operator>)
    {
        // sum: accumulate from all overlapping tile ghost boxes (resize already zero-inits)
        for (auto const& tile : src.field())
        {
            auto const tile_gb = tile.layout().AMRGhostBoxFor(pq);
            if (auto const overlap = amr_selection_box * tile_gb)
                for (auto const& bix : *overlap)
                {
                    auto const dst_lix = (bix - amr_selection_box.lower).as_unsigned();
                    auto const src_lix = (bix - tile_gb.lower).as_unsigned();
                    Operator{view(dst_lix)}(tile()(src_lix));
                }
        }
    }
    else
    {
        // set/max: use reduce_single — one canonical tile per cell, correct for all values
        reduce_single(view, src.field, amr_selection_box);
    }
}

template<typename Operator, typename... Args>
void append_from_fields(FieldBox<GridTileSet<Args...> const> const& src, auto& vec)
{
    append_from_fields<Operator>(
        FieldBox<FieldTileSet<Args...> const>{src.field, src.amr_box, src.lcl_box}, vec);
}

template<typename Operator, typename... Args>
void append_from_fields(FieldBox<Args...> const& src, auto& vec)
{
    // Operator unused assumed SetEqual
    // reserve vec before use!
    auto src_it = src.lcl_box.begin();
    for (; src_it != src.lcl_box.end(); ++src_it)
        vec.push_back(src.field(*src_it));
}

template<typename Field_t>
template<typename Operator>
void FieldBox<Field_t>::append_to(std::vector<value_type>& vec) const
{
    append_from_fields<Operator>(*this, vec);
}

template<typename... T0s, typename... T1s>
void copy_fields(FieldTileSet<T0s...>& dst, FieldTileSet<T1s...> const& src)
{
    for (std::size_t idx = 0; idx < src().size(); ++idx)
        copy_fields(dst[idx](), src[idx]());
}

template<typename... T0s, auto opts>
void copy_fields(FieldTileSet<T0s...>& dst, basic::Field<opts> const& src)
{
    PHARE_LOG_SCOPE(3, "copy_fields<FieldTileSet,basic::Field>");

    assert(dst().size());
    auto const& patch_layout = dst()[0].layout().copy_as(dst.box());

    for (auto& tile : dst)
        FieldBox{tile(), tile.layout(), tile.ghost_box()}.op(
            FieldBox{src, patch_layout, tile.ghost_box()});
}

template<typename... T0s, auto opts>
void copy_fields(basic::Field<opts>& dst, FieldTileSet<T0s...> const& src)
{
    reduce_single_<SetEqual<typename basic::Field<opts>::value_type>>(src, dst);
}

template<auto opts0, auto opts1>
void copy_fields(basic::Field<opts0>& dst, basic::Field<opts1> const& src)
{
    std::memcpy(dst.data(), src.data(),
                src.size() * sizeof(typename basic::Field<opts0>::value_type));
}

template<typename Operator, typename Grid_t, typename GridTiles_t>
auto& reduce_single_(GridTiles_t const& tiles, Grid_t& grid)
    requires(is_field_tile_set_v<GridTiles_t> && is_field_v<Grid_t>)
{
    auto const pq           = tiles.physicalQuantity();
    auto const patch_layout = tiles[0].layout().copy_as(tiles.box());
    auto const patch_box    = patch_layout.AMRBox();
    auto const ghost_box    = patch_layout.AMRGhostBoxFor(pq);
    auto const expander     = BoxExpander{ghost_box, patch_box};

    grid.reshape(tiles.shape());
    grid.zero();

    for (auto const& tile : tiles())
    {
        auto const tile_gb              = tile.layout().AMRGhostBoxFor(pq);
        auto const is_border            = *(tile_gb * patch_box) != tile_gb;
        decltype(tile_gb) const& tile_b = *tile;
        for (auto const& bix : is_border ? expander(tile_b) : tile_b)
        {
            auto const lix      = (bix - ghost_box.lower).as_unsigned();
            auto const tile_lix = (bix - tile_gb.lower).as_unsigned();
            grid(lix)           = tile()(tile_lix);
        }
    }
    return grid;
}

template<typename Operator, typename Grid_t, typename GridTiles_t>
auto& reduce_single_(Grid_t const& grid, GridTiles_t& tiles)
    requires(is_field_v<Grid_t> && is_field_tile_set_v<GridTiles_t>)
{
    auto const pq           = tiles.physicalQuantity();
    auto const get_box      = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    auto const patch_layout = tiles[0].layout().copy_as(tiles.box());
    auto const patch_box    = patch_layout.AMRBox();
    auto const ghost_box    = patch_layout.AMRGhostBoxFor(pq);
    auto const expander     = BoxExpander{ghost_box, patch_box};
    tiles.zero();

    assert(array_equals(tiles.shape(), grid.shape()));

    for (auto& tile : tiles())
    {
        auto const tile_gb              = tile.layout().AMRGhostBoxFor(pq);
        auto const is_border            = *(tile_gb * patch_box) != tile_gb;
        decltype(tile_gb) const& tile_b = *tile;
        for (auto const& bix : is_border ? expander(tile_b) : tile_b)
        {
            auto const lix      = (bix - ghost_box.lower).as_unsigned();
            auto const tile_lix = (bix - tile_gb.lower).as_unsigned();
            tile()(tile_lix)    = grid(lix);
        }
    }
    return grid;
}

template<typename Dst, typename Src>
auto& reduce_single(Dst& dst, Src const& src)
{
    if constexpr (is_field_tile_set_v<Src>)
        reduce_single_<SetEqual<typename Dst::value_type>>(src, dst);
    else if constexpr (is_field_tile_set_v<Dst>)
        reduce_single_<SetEqual<typename Dst::value_type>>(src, dst);
    else
    {
        static_assert(is_field_v<Src> and is_field_v<Dst>);
        copy_fields(dst, src);
    }
    return dst;
}

template<typename Tiles>
auto reduce_single(Tiles const& input)
    requires(is_field_tile_set_v<Tiles>)
{
    using Grid_t = Tiles::grid_type;
    Grid_t grid{input.name(), input.physicalQuantity(), input.shape()};
    reduce_single(grid, input);
    return grid;
}

template<typename Tiles>
auto& reduce_single(Tiles const& input)
    requires(!is_field_tile_set_v<Tiles>)
{
    return input;
}

template<typename Operator, typename Grid_t, typename GridTiles_t>
auto& reduce_into_(GridTiles_t const& tiles, Grid_t& grid)
{
    PHARE_LOG_SCOPE(3, "reduce_into<GridTileSet,Grid>");

    grid.reshape(tiles.shape());
    grid.zero();
    auto const& patch_layout = tiles[0].layout().copy_as(tiles.box());

    auto const pq      = tiles.physicalQuantity();
    auto const get_box = [&](auto const& tile) { return tile.layout().AMRGhostBoxFor(pq); };
    for (auto const& tile : tiles())
    {
        auto const& tile_layout = tile.layout();
        auto const& tile_box    = get_box(tile);
        FieldBox{grid, patch_layout, patch_layout.AMRToLocal(tile_box)}. //
            template op<Operator>(core::FieldBox{tile(), tile_layout, tile_box});
    }
    return grid;
}

template<typename Dst, typename Src>
auto& reduce_into(Dst& dst, Src const& src)
{
    if constexpr (is_field_tile_set_v<Src>)
        reduce_into_<PlusEquals<typename Dst::value_type>>(src, dst);
    else
    {
        static_assert(is_field_v<Dst> and is_field_v<Src>);
        copy_fields(dst, src);
    }
    return dst;
}

template<typename Dst, typename TiledField>
auto& reduce_single(Dst& dst, TiledField const& tiles,
                    Box<int, TiledField::dimension> const& selection)
    requires(is_field_tile_set_v<TiledField>)
{
    if (tiles().size() == 0)
        return dst;

    auto const pq          = tiles.physicalQuantity();
    auto const& first_tile = tiles()[0];
    auto const patch_gb    = first_tile.layout().copy_as(tiles.box()).AMRGhostBoxFor(pq);
    auto const patch_fb    = shrink(patch_gb, first_tile.layout().options.field_ghost_width);
    BoxExpander const expander{patch_gb, patch_fb};

    for (auto const& tile : tiles())
    {
        auto const tile_gb     = tile.layout().AMRGhostBoxFor(pq);
        auto const expanded_fb = expander(tile.field_box());
        if (auto const overlap = selection * expanded_fb)
        {
            for (auto const& bix : *overlap)
            {
                auto const dst_lix = (bix - selection.lower).as_unsigned();
                auto const src_lix = (bix - tile_gb.lower).as_unsigned();
                dst(dst_lix)       = tile()(src_lix);
            }
        }
    }
    return dst;
}

template<typename Tiles>
auto reduce(Tiles const& input)
    requires(is_field_tile_set_v<Tiles>)
{
    using Grid_t = Tiles::grid_type;
    Grid_t grid{input.name(), input.physicalQuantity(), input.shape()};
    reduce_into(grid, input);
    return grid;
}

template<typename Tiles>
auto& reduce(Tiles const& input)
    requires(!is_field_tile_set_v<Tiles>)
{
    return input;
}


template<typename Field_t>
void fill_domain_box(Field_t& field, auto const& layout, auto const v)
    requires(!is_field_tile_set_v<Field_t>)
{
    auto const pq = field.physicalQuantity();
    for (auto const& bix : layout.domainBoxFor(pq))
        field(bix) = v;
}

template<typename Field_t>
void fill_domain_box(Field_t& field, auto const& layout, auto const v)
    requires(is_field_tile_set_v<Field_t>)
{
    auto const pq = field.physicalQuantity();
    for (auto& tile : field())
        for (auto const& bix : tile.layout().domainBoxFor(pq))
            tile()(bix) = v;
}

template<typename Field_t>
void fill_ghost(Field_t& field, auto const& layout, auto const v)
    requires(!is_field_tile_set_v<Field_t>)
{
    field[NdArrayMask{0, layout.options.field_ghost_width}] = v;
}

template<typename Field_t>
void fill_ghost(Field_t& field, auto const& layout, auto const v)
    requires(is_field_tile_set_v<Field_t>)
{
    auto const pq        = field.physicalQuantity();
    auto const patch_box = layout.AMRBoxFor(field);
    for (auto& tile : field())
    {
        auto const tile_gb = tile.layout().AMRGhostBoxFor(pq);
        if (auto const is_border = *(tile_gb * patch_box) != tile_gb; not is_border)
            continue;

        for (auto const& bix : tile_gb)
            if (!isIn(bix, patch_box))
                tile()((bix - tile_gb.lower).as_unsigned()) = v;
    }
}

} // namespace PHARE::core

#endif
