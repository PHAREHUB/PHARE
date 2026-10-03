#ifndef PHARE_CORE_DATA_FIELD_FIELD_TILES_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_TILES_HPP

#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"

#include <tuple>
#include <string>
#include <cstddef>
#include <stdexcept>
#include <algorithm>
#include <type_traits>

namespace PHARE::core::basic
{
template<typename Field_t, std::size_t rank>
struct TensorField;
}

namespace PHARE::core
{
template<typename GridLayout_t, typename Grid_t, typename Field_t>
struct GridTile;

template<typename GridLayout_t, typename Field_t>
class FieldTile : public Box<std::int32_t, GridLayout_t::dimension>
{
public:
    auto constexpr static dimension  = GridLayout_t::dimension;
    auto constexpr static alloc_mode = Field_t::alloc_mode;
    using Super                      = Box<std::int32_t, dimension>;
    using value_type                 = Field_t;
    using type                       = Field_t::type;

    template<typename... Args>
    FieldTile(GridLayout_t const& layout, Field_t const& field)
        : Super{layout.AMRBox()}
        , field_{field}
        , layout_{layout}
    {
    }

    auto& operator()() { return field_; }
    auto& operator()() const { return field_; }
    Super& operator*() { return *this; }
    Super const& operator*() const { return *this; }
    auto& layout() const { return layout_; }

    auto ghost_box() const { return ghost_box_; }
    auto field_box() const { return shrink(ghost_box(), GridLayout_t::options.field_ghost_width); }

    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point)
    {
        return field_(point);
    }
    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) const
    {
        return field_(point);
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == dimension)

    {
        return field_(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == dimension)

    {
        return field_(to_point<std::uint32_t>(ijk...));
    }

    NO_DISCARD auto physicalQuantity() const { return field_.physicalQuantity(); }

    bool isUsable() const { return field_.isUsable(); }
    bool isSettable() const { return !isUsable(); }

    auto& links() { return _links; }
    auto& links() const { return _links; }
    auto& link(std::size_t const idx) { return _links[idx]; }

private:
    Field_t field_;
    GridLayout_t layout_;
    Box<int, dimension> ghost_box_   = layout_.AMRGhostBoxFor(physicalQuantity());
    std::array<FieldTile*, 7> _links = ConstArray<FieldTile*, 7>(nullptr);
};

template<typename GridLayout_t, typename Grid_t, typename Field_t>
struct FieldTileSetter
{
    using Tile_t     = GridTile<GridLayout_t, Grid_t, Field_t>;
    using Tile_vt    = FieldTile<GridLayout_t, Field_t>;
    using Span_t     = ViewSpan<Tile_vt, Tile_t>;
    using Span_pt    = ViewSpan<Tile_vt*, Tile_t*>;
    using value_type = TileSetView<Tile_vt, Span_t, Span_pt>;
};

} // namespace PHARE::core

namespace PHARE::core::basic
{
template<typename GridLayout_, typename Grid_t, typename Field_t>
class FieldTileSet : public FieldTileSetter<GridLayout_, Grid_t, Field_t>::value_type
{
    using local_types = FieldTileSetter<GridLayout_, Grid_t, Field_t>;

public:
    using GridLayout_t               = GridLayout_;
    using grid_type                  = Grid_t;
    using field_type                 = Field_t;
    using Super                      = local_types::value_type;
    using value_type                 = local_types::Tile_vt;
    using physical_quantity_type     = Grid_t::physical_quantity_type;
    using type                       = Grid_t::type;
    auto constexpr static dimension  = GridLayout_t::dimension;
    auto constexpr static alloc_mode = Grid_t::alloc_mode;

    FieldTileSet(auto physicalQuantity)
        : Super{{}, nullptr, 0, nullptr, {}}
        , qty_{physicalQuantity}
    {
    }

    FieldTileSet(FieldTileSet const&)            = default;
    FieldTileSet(FieldTileSet&&)                 = delete;
    FieldTileSet& operator=(FieldTileSet const&) = delete;
    FieldTileSet& operator=(FieldTileSet&&)      = delete;

    void setBuffer(std::nullptr_t ptr) { super() = Super{{}, nullptr, 0, nullptr, {}}; }

    template<typename FieldLike>
    void setBuffer(FieldLike* const field)
    {
        auto data = field ? field->data() : nullptr;
        if (data)
        {
            super()        = field->super().as([](auto&&... args) mutable {
                auto&& [box, tiles_data, n_tiles, cells_data, cells_shape]
                    = std::forward_as_tuple(args...);
                assert(cells_data);
                return Super{box, &*tiles_data[0], n_tiles,
                             reinterpret_cast<FieldTile<GridLayout_t, Field_t>**>(cells_data),
                             cells_shape};
            });
            ghost_box_     = field->layout().AMRGhostBoxFor(qty_);
            max_tile_size_ = field->max_tile_size();
        }
        else
            super() = Super{{}, nullptr, 0, nullptr, {}};
    }

    void zero()
    {
        if (!Super::data())
            throw std::runtime_error("invalid state");

        for (auto& tile : *this)
            std::fill(tile().data(), tile().data() + tile().size(), 0);
    }

    void fill(auto const v)
    {
        if (!Super::data())
            throw std::runtime_error("invalid state");
        for (auto& tile : *this)
            std::fill(tile().data(), tile().data() + tile().size(), v);
    }

    auto& operator()() { return super()(); }
    auto& operator()() const { return super()(); }

    bool isUsable() const { return Super::data() != nullptr; }
    bool isSettable() const { return !isUsable(); }

    void copyData(FieldTileSet const& that)
    {
        for (std::size_t tidx = 0; tidx < super().size(); ++tidx)
            std::copy(that[tidx]().data(), that[tidx]().data() + that[tidx]().size(),
                      super()[tidx]().data());
    }
    NO_DISCARD auto size() const { return ghost_box_.size(); } // NOT ntiles!
    NO_DISCARD auto shape() const { return *ghost_box_.shape().as_unsigned(); }
    NO_DISCARD auto ghost_box() const { return ghost_box_; }

    NO_DISCARD auto at(auto const&... args) { return super().at(args...); }
    NO_DISCARD auto at(auto const&... args) const { return super().at(args...); }

    void sync_inner_ghosts()
    {
        visit_tile_neighbours(super(), [](auto& t0, auto& t1) {
            if (auto const overlap = t0.ghost_box() * t1.field_box())
                for (auto const& bix : *overlap)
                {
                    auto const& t0_lix = (bix - t0.ghost_box().lower).as_unsigned();
                    auto const& t1_lix = (bix - t1.ghost_box().lower).as_unsigned();
                    t0()(t0_lix)       = t1()(t1_lix);
                }
        });
    }


    void notZero() const
    {
        for (auto& tile : super())
            for (auto const& e : tile())
                if (std::abs(e) < 1e-15)
                    throw std::runtime_error("ZERO");
    }


    auto max_tile_size() const { return max_tile_size_; }
    auto ntiles() const { return Super::size(); }

    NO_DISCARD auto physicalQuantity() const { return qty_; }

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    physical_quantity_type qty_;
    Box<int, dimension> ghost_box_{};
    std::uint32_t max_tile_size_;
};


} // namespace PHARE::core::basic

namespace PHARE::core
{

template<typename GridLayout_t, typename Grid_t, typename Field_t>
class FieldTileSet : public basic::FieldTileSet<GridLayout_t, Grid_t, Field_t>
{
public:
    using Super = basic::FieldTileSet<GridLayout_t, Grid_t, Field_t>;

    FieldTileSet(std::string const& name, auto physicalQuantity)
        : Super{physicalQuantity}
        , name_{name}
    {
    }

    NO_DISCARD auto& name() const { return name_; }

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    std::string const name_;
};


template<typename T>
concept is_field_tile_set_c = requires(T* p) {
    []<typename GL, typename Arr, typename PQ>(basic::FieldTileSet<GL, Arr, PQ> const*) {}(p);
};

template<typename T>
inline constexpr bool is_field_tile_set_v = is_field_tile_set_c<T>;

template<typename T>
concept has_tiled_field_type_c = requires {
    typename std::remove_cvref_t<T>::field_type;
} and is_field_tile_set_v<typename std::remove_cvref_t<T>::field_type>;

template<typename T>
auto tile_at(T&& tiled, std::size_t const i)
{
    using Tiled_t = std::remove_cvref_t<T>;
    if constexpr (is_field_tile_set_v<Tiled_t>)
        return tiled()[i]();
    else if constexpr (has_tiled_field_type_c<T>)
    {
        using Field_t = typename Tiled_t::field_type;
        using Tile_t  = typename Field_t::value_type::value_type;
        using V_t     = basic::TensorField<Tile_t, Tiled_t::rank>;
        return tiled.template as<V_t>([&](auto& c) { return c()[i](); });
    }
    else
        return std::forward<T>(tiled);
}

template<typename... Tiled>
auto tiles_at(std::size_t const i, Tiled&&... tiled)
{
    return std::make_tuple(tile_at(std::forward<Tiled>(tiled), i)...);
}

template<typename T>
auto tile_count(T const& tiled)
{
    if constexpr (is_field_tile_set_v<std::remove_cvref_t<T>>)
        return tiled().size();
    else
        return tiled[0]().size();
}

template<typename T>
auto& tile_layout(T const& tiled, std::size_t const i)
{
    if constexpr (is_field_tile_set_v<std::remove_cvref_t<T>>)
        return tiled()[i].layout();
    else // tensorfield expected
        return tiled[0]()[i].layout();
}

template<typename Fn, typename Tiled0, typename... Args>
void tile_exec_with_layout(Fn&& fn, Tiled0&& tiled0, Args&&... args)
{
    for (std::size_t i = 0; i < tile_count(tiled0); ++i)
        fn(tile_layout(tiled0, i), tile_at(tiled0, i), tile_at(args, i)...);
}

void sync_inner_ghosts(auto& vf)
{
    if constexpr (is_field_tile_set_v<std::remove_cvref_t<decltype(vf)>>)
        vf.sync_inner_ghosts();
    else
        for (auto& c : vf)
            c.sync_inner_ghosts();
}


} // namespace PHARE::core

#endif // PHARE_CORE_DATA_FIELD_FIELD_TILES_HPP
