#ifndef PHARE_CORE_DATA_FIELD_FIELD_OVERLAPS_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_OVERLAPS_HPP


// #include "core/data/grid/gridlayout.hpp"
#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"
// #include "core/data/field/field_box.hpp"
// #include "core/data/grid/grid_tiles.hpp"
#include "core/utilities/box/box_span.hpp"
#include "core/data/tiles/tile_set_overlaps.hpp"


#include <map>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <cstring>



namespace PHARE::core
{


template<typename GridLayout_, auto opts>
struct FieldTileOverlaps
{
    using Field_t     = basic::Field<opts>;
    using FieldTile_t = FieldTile<GridLayout_, Field_t>;
    using Quantity    = decltype(opts)::physical_quantity_type;


    // when you're sure it should exist
    static auto& getQuantity(auto const& layout, Quantity const& pq)
    {
        return levels.at(layout.levelNumber()).patches.at(to_string(layout.AMRBox())).at(pq);
    }
    // when you're sure it should exist
    template<typename Tiles>
    static auto& getQuantity(auto const& layout, Tiles const& tiles)
        requires(has_physicalQuantity_v<Tiles>)
    {
        return getQuantity(layout, tiles.physicalQuantity());
    }

    // when you're not sure it should exist
    static auto& getOrCreateQuantity(auto const& layout, auto& field);


    static auto build_inner_tile_overlaps(auto& field);
    static void sync_inner_ghosts(auto& field, auto const& overlaps);
    static void reset(int lvl) { levels[lvl].clear(); }
    static void reset() { levels.clear(); }


    struct Overlap
    {
        Overlap(auto _src, auto const& d, auto const& s)
            : src{_src}
            , lcl_dst{d}
            , lcl_src{s}
        {
        }

        Field_t const* const src;
        Box<std::uint32_t, opts.dimension> lcl_dst, lcl_src;
    };

    struct Level
    {
        struct Patch
        {
            using TileSpan_t = NdSpanSet<opts.dimension, FieldTile_t*>;

            struct Tile
            {
                Tile(auto& field_tile)
                    : domain_slabs{make_box_span(field_tile.layout().domainBoxFor(field_tile()))}
                    , ghost_slabs{make_box_span(field_tile.layout().ghostBoxFor(field_tile()))}
                {
                }

                BoxSpan<opts.dimension> domain_slabs, ghost_slabs;
                std::vector<Overlap> overlaps{};
            };

            Patch(auto const& layout, auto& field)
                : tile_span{make_qty_nd_span_set_from(field, layout, field.physicalQuantity())}
            {
                static_assert(not std::is_const_v<decltype(field)>);
            }

            NO_DISCARD auto& operator[](std::size_t const i) { return tiles[i]; }
            NO_DISCARD auto& operator[](std::size_t const i) const { return tiles[i]; }
            NO_DISCARD auto begin() { return tiles.begin(); }
            NO_DISCARD auto begin() const { return tiles.begin(); }
            NO_DISCARD auto end() { return tiles.end(); }
            NO_DISCARD auto end() const { return tiles.end(); }

            TileSpan_t tile_span;
            std::vector<Tile> tiles{};
        };

        std::map<std::string, std::map<Quantity, Patch>> patches{};
    };

    static inline std::map<int, Level> levels{};
};


template<typename GL, auto opts>
auto& FieldTileOverlaps<GL, opts>::getOrCreateQuantity(auto const& layout, auto& field)
{
    static_assert(not std::is_const_v<decltype(field)>);

    using Patch_t = Level::Patch;

    if (!levels.count(layout.levelNumber()))
        levels.try_emplace(layout.levelNumber());

    auto& level           = levels.at(layout.levelNumber());
    auto const& patch_key = to_string(layout.AMRBox());
    if (!level.patches.count(patch_key))
        level.patches.try_emplace(patch_key);

    auto& patch_map = level.patches.at(patch_key);
    if (!patch_map.count(field.physicalQuantity()))
        patch_map.try_emplace(field.physicalQuantity(), Patch_t{layout, field});

    auto& patch = patch_map.at(field.physicalQuantity());
    if (!patch.tiles.size())
        patch.tiles = build_inner_tile_overlaps(field);

    return patch;
}

template<typename GL, auto opts>
auto FieldTileOverlaps<GL, opts>::build_inner_tile_overlaps(auto& field)
{
    using TileOverlap = Level::Patch::Tile;

    std::vector<TileOverlap> tiles;
    tiles.reserve(field().size());

    for (auto& tile : field())
        tiles.emplace_back(tile);

    for (std::size_t ti0 = 0; ti0 < field().size() - 1; ++ti0)
        for (std::size_t ti1 = ti0 + 1; ti1 < field().size(); ++ti1)
        {
            auto& t0 = field()[ti0];
            auto& t1 = field()[ti1];

            if (auto const t0_overlap = t0.ghost_box() * t1.field_box())
                tiles[ti0].overlaps.emplace_back(
                    &t1(), as_unsigned(shift(*t0_overlap, t0.ghost_box().lower * -1)),
                    as_unsigned(shift(*t0_overlap, t1.ghost_box().lower * -1)));
            else
                continue; // second can't overlap if first doesn't

            if (auto const t1_overlap = t1.ghost_box() * t0.field_box())
                tiles[ti1].overlaps.emplace_back(
                    &t0(), as_unsigned(shift(*t1_overlap, t1.ghost_box().lower * -1)),
                    as_unsigned(shift(*t1_overlap, t0.ghost_box().lower * -1)));
        }


    return tiles;
}



} // namespace PHARE::core


#endif
