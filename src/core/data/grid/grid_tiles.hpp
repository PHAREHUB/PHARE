#ifndef PHARE_CORE_DATA_GRID_GRID_TILES_HPP
#define PHARE_CORE_DATA_GRID_GRID_TILES_HPP

#include "core/def/phare_config.hpp"

#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/data/field/field.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/tiles/tile_set.hpp"
#include "core/data/tiles/tile_set_traversal.hpp"


#include <optional>
#include <tuple>
#include <string>
#include <cstddef>
#include <stdexcept>


namespace PHARE::core
{
template<typename Field_t>
class FieldBox;

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

    auto& operator()() _PHARE_ALL_FN_ { return field_; }
    auto& operator()() const _PHARE_ALL_FN_ { return field_; }
    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }
    auto& layout() const _PHARE_ALL_FN_ { return layout_; }

    auto ghost_box() const _PHARE_ALL_FN_ { return ghost_box_; }
    auto field_box() const _PHARE_ALL_FN_ { return shrink(ghost_box(), GridLayout_t::nbrGhosts()); }


    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) _PHARE_ALL_FN_
    {
        return field_(point);
    }
    template<template<typename, std::size_t> typename Point_t>
    auto& operator()(Point_t<std::uint32_t, dimension> const& point) const _PHARE_ALL_FN_
    {
        return field_(point);
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk)
        requires(sizeof...(IJK) == dimension)
    _PHARE_ALL_FN_
    {
        return field_(to_point<std::uint32_t>(ijk...));
    }
    template<typename... IJK>
    auto& operator()(IJK const&... ijk) const
        requires(sizeof...(IJK) == dimension)
    _PHARE_ALL_FN_
    {
        return field_(to_point<std::uint32_t>(ijk...));
    }

    NO_DISCARD auto physicalQuantity() const _PHARE_ALL_FN_ { return field_.physicalQuantity(); }

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



template<bool hasGhosts = true>
auto static grid_cells(auto&&... args)
{
    auto const& [layout, qty] = std::forward_as_tuple(args...);

    if constexpr (hasGhosts)
        return layout.allocSize(qty);
    else
        return *(layout.AMRBox().shape().as_unsigned() - 1);
}



template<typename GridLayout_t, typename Grid_t, typename Field_t>
struct GridTile : public FieldTile<GridLayout_t, Field_t>
{
    auto constexpr static dimension = GridLayout_t::dimension;
    using Super                     = FieldTile<GridLayout_t, Field_t>;
    using value_type                = Grid_t::Super;
    using NdArray_t                 = Grid_t::Super;
    using Super::operator();

    GridTile(auto const& layout, auto const& pq)
        : Super{layout, Field_t{pq}}
        , arr{grid_cells(layout, pq)}
    {
        reset();
    }

    GridTile(auto const& box, auto const& patch_layout,
             auto const& pq) // called when making tiles
        : GridTile(patch_layout.copy_as(box), pq)
    {
    }

    GridTile(GridTile const& that)
        : Super{that.layout(), Field_t{(*that).physicalQuantity()}}
        , arr{that.arr}
    {
        reset();
    }

    GridTile(GridTile&&)                 = default;
    GridTile& operator=(GridTile const&) = delete;
    GridTile& operator=(GridTile&&)      = default;

    Super& operator*() _PHARE_ALL_FN_ { return *this; }
    Super const& operator*() const _PHARE_ALL_FN_ { return *this; }

    NO_DISCARD auto physicalQuantity() const _PHARE_ALL_FN_ { return (**this).physicalQuantity(); }

private:
    void reset() { (**this)() = Field_t{(**this).physicalQuantity(), arr.data(), arr.shape()}; }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

    NdArray_t arr;
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

    FieldTileSet(auto physicalQuantity) _PHARE_ALL_FN_ : Super{{}, nullptr, 0, nullptr, {}},
                                                         qty_{physicalQuantity}
    {
    }


    FieldTileSet(FieldTileSet const&) _PHARE_ALL_FN_            = default;
    FieldTileSet(FieldTileSet&&) _PHARE_ALL_FN_                 = delete;
    FieldTileSet& operator=(FieldTileSet const&) _PHARE_ALL_FN_ = delete;
    FieldTileSet& operator=(FieldTileSet&&) _PHARE_ALL_FN_      = delete;


    void setBuffer(std::nullptr_t ptr) _PHARE_ALL_FN_
    {
        super() = Super{{}, nullptr, 0, nullptr, {}};
    }

    template<typename FieldLike>
    void setBuffer(FieldLike* const field) _PHARE_ALL_FN_
    {
        auto data = field ? field->data() : nullptr;
        if (data)
        {
            super() = field->super().as([](auto&&... args) mutable {
                auto&& [box, tiles_data, n_tiles, cells_data, cells_shape]
                    = std::forward_as_tuple(args...);
                assert(cells_data);
                return Super{box, &*tiles_data[0], n_tiles,
                             reinterpret_cast<FieldTile<GridLayout_t, Field_t>**>(cells_data),
                             cells_shape};
            });
            // check();
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

        using vec_helper = PHARE::Vector<type, alloc_mode>;
        for (auto& tile : *this)
            vec_helper::fill(tile().data(), tile().size(), 1e-36);
    }


    void fill(auto const v)
    {
        if (!Super::data())
            throw std::runtime_error("invalid state");

        using vec_helper = PHARE::Vector<type, alloc_mode>;
        for (auto& tile : *this)
            vec_helper::fill(tile().data(), tile().size(), v);
    }


    auto& operator()() _PHARE_ALL_FN_
    {
        PHARE_ASSERT(isUsable());
        return super()();
    }
    auto& operator()() const _PHARE_ALL_FN_
    {
        PHARE_ASSERT(isUsable());
        return super()();
    }


    void setData(auto* const data) _PHARE_ALL_FN_
    {
        throw std::runtime_error("fix?");
        /*Super::setBuffer(data);*/
    }

    bool isUsable() const
    {
        auto const b = Super::data() != nullptr;
        // if (b)
        //     check();
        return b;
    }
    bool isSettable() const { return !isUsable(); }

    // auto& operator*() { return *this; }       // interop atm
    // auto& operator*() const { return *this; } // interop atm

    void copyData(FieldTileSet const& that)
    {
        for (std::size_t tidx = 0; tidx < super().size(); ++tidx)
            std::copy(that[tidx]().data(), that[tidx]().data() + that[tidx]().size(),
                      super()[tidx]().data());
    }
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return ghost_box_.size(); } // NOT ntiles!
    NO_DISCARD auto shape() const _PHARE_ALL_FN_ { return *ghost_box_.shape().as_unsigned(); }
    NO_DISCARD auto ghost_box() const _PHARE_ALL_FN_ { return ghost_box_; }


    NO_DISCARD auto at(auto const&... args) { return super().at(args...); }
    NO_DISCARD auto at(auto const&... args) const { return super().at(args...); }


    void check(bool const nan = false) const
    {
        // for (auto& tile : super())
        // {
        //     assert(tile().size() < static_cast<std::size_t>(1e6));
        //     tile().check();
        // }
    }

    void sync_inner_ghosts()
    {
        // KUL_DBG_FUNC_ENTER
        // PHARE_LOG_LINE_SS("");
        // visit_tile_neighbours(super(), [](auto& t0, auto& t1) {
        //     //
        //     if (auto const t0_overlap = t0.ghost_box() * t1.field_box())
        //         for (auto const& bix : *t0_overlap)
        //         {
        //             auto const& t0_lix = (bix - t0.ghost_box().lower).as_unsigned();
        //             auto const& t1_lix = (bix - t1.ghost_box().lower).as_unsigned();
        //             t0()(t0_lix)       = t1()(t1_lix);
        //         }

        //     // if (auto const t1_overlap = t1.ghost_box() * t0.field_box())
        //     //     for (auto const& bix : *t1_overlap)
        //     //     {
        //     //         auto const& t0_lix = (bix - t0.ghost_box().lower).as_unsigned();
        //     //         auto const& t1_lix = (bix - t1.ghost_box().lower).as_unsigned();
        //     //         t1()(t1_lix)       = t0()(t0_lix);
        //     //     }
        // });

        for (std::size_t ti0 = 0; ti0 < super().size() - 1; ++ti0)
        {
            for (std::size_t ti1 = ti0 + 1; ti1 < super().size(); ++ti1)
            {
                auto& t0 = super()[ti0];
                auto& t1 = super()[ti1];

                if (auto const t0_overlap = t0.ghost_box() * t1.field_box())
                    for (auto const& bix : *t0_overlap)
                    {
                        auto const& t0_lix = (bix - t0.ghost_box().lower).as_unsigned();
                        auto const& t1_lix = (bix - t1.ghost_box().lower).as_unsigned();
                        t0()(t0_lix)       = t1()(t1_lix);
                    }

                if (auto const t1_overlap = t1.ghost_box() * t0.field_box())
                    for (auto const& bix : *t1_overlap)
                    {
                        auto const& t0_lix = (bix - t0.ghost_box().lower).as_unsigned();
                        auto const& t1_lix = (bix - t1.ghost_box().lower).as_unsigned();
                        t1()(t1_lix)       = t0()(t0_lix);
                    }
            }
        }
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

    template<typename, typename, typename>
    friend std::ostream& operator<<(std::ostream&, FieldTileSet const&);

    NO_DISCARD auto physicalQuantity() const _PHARE_ALL_FN_ { return qty_; }

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
    // using grid_type  = Grid_t;
    // using field_type = Field_t;
    using Super = basic::FieldTileSet<GridLayout_t, Grid_t, Field_t>;
    // using value_type                 = local_types::Tile_vt;
    // using physical_quantity_type     = Grid_t::physical_quantity_type;
    // using type                       = Grid_t::type;
    // auto constexpr static dimension  = GridLayout_t::dimension;
    // auto constexpr static alloc_mode = Grid_t::alloc_mode;

    FieldTileSet(std::string const& name, auto physicalQuantity) _PHARE_ALL_FN_
        : Super{physicalQuantity},
          name_{name}
    {
    }

    // FieldTileSet(FieldTileSet const&) _PHARE_ALL_FN_            = default;
    // FieldTileSet(FieldTileSet&&) _PHARE_ALL_FN_                 = delete;
    // FieldTileSet& operator=(FieldTileSet const&) _PHARE_ALL_FN_ = delete;
    // FieldTileSet& operator=(FieldTileSet&&) _PHARE_ALL_FN_      = delete;


    NO_DISCARD auto& name() const { return name_; }

private:
    Super& super() { return *this; }
    Super const& super() const { return *this; }

    std::string const name_;
};


template<typename GridLayout_t, typename Grid_t, typename Field_t>
class GridTileSet : public TileSet<GridTile<GridLayout_t, Grid_t, Field_t>, Grid_t::alloc_mode>,
                    public FieldTileSet<GridLayout_t, Grid_t, Field_t>

{
public:
    using grid_type  = Grid_t;
    using field_type = Field_t;
    using Super      = TileSet<GridTile<GridLayout_t, Grid_t, Field_t>, Grid_t::alloc_mode>;
    using View       = FieldTileSet<GridLayout_t, Grid_t, Field_t>;
    using type       = Grid_t::value_type;
    using physical_quantity_type = Grid_t::physical_quantity_type;
    using value_type             = GridTile<GridLayout_t, Grid_t, Field_t>;
    using Super::operator[];
    using View::operator();
    using View::box;

    auto constexpr static alloc_mode = Grid_t::alloc_mode;
    auto constexpr static dimension  = Grid_t::dimension;

    GridTileSet(std::string const& name, GridLayout_t const& layout,
                physical_quantity_type const qty, std::optional<type> val = std::nullopt)
        : Super{layout.AMRBox(), layout, qty}
        , View{name, qty}
        , name_{name}
        , layout_{layout}
        , max_tile_size_{get_max_tile_size()}
    {
        Super::build_links(super());
        assert(this->super()[0]().data());
        View::setBuffer(this);
        assert(View::isUsable());
        assert((**this)[0]().data());
        if (val)
            fill(*val);
        // ptr = &TileOverlaps_t::getOrCreateQuantity(layout_, *this);
    }

    GridTileSet(GridTileSet const& that)
        : Super{that.super()}
        , View{that}
        , layout_{that.layout_}
        , max_tile_size_{get_max_tile_size()}
    {
        Super::build_links(super());
        View::setBuffer(this);
        assert((**this)[0]().data());
        assert(View::isUsable());
        // ptr = &TileOverlaps_t::getOrCreateQuantity(layout_, *this);
    }

    GridTileSet(GridTileSet&&) = default;

    GridTileSet& operator=(GridTileSet const&) = delete;
    GridTileSet& operator=(GridTileSet&&)      = delete;


    template<typename T>
    void fill(T const v)
    {
        for (auto& tile : *this)
            tile().fill(v);
    }


    auto& operator()() { return super()(); }
    auto& operator()() const { return super()(); }

    NO_DISCARD auto at(auto const&... args) { return super().at(args...); }
    NO_DISCARD auto at(auto const&... args) const { return super().at(args...); }

    View& operator*() { return *this; }
    View const& operator*() const { return *this; }

    auto data() { return super().data(); }
    auto data() const { return super().data(); }
    NO_DISCARD auto begin() { return super().begin(); }
    NO_DISCARD auto begin() const { return super().begin(); }
    NO_DISCARD auto end() { return super().end(); }
    NO_DISCARD auto end() const { return super().end(); }
    NO_DISCARD auto size() const _PHARE_ALL_FN_ { return View::size(); }
    auto& layout() const { return layout_; }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

    NO_DISCARD auto& name() const { return name_; }
    auto max_tile_size() const { return max_tile_size_; }

private:
    auto get_max_tile_size() const
    {
        std::uint32_t size = 0;
        for (auto const& tile : super()())
            if (std::uint32_t tile_size = product(tile.ghost_box().shape()); tile_size > size)
                size = tile_size;
        return size;
    }

    std::string name_;
    GridLayout_t layout_;
    std::uint32_t max_tile_size_;

    // FieldTileOverlaps<field_opts>::Level::Patch* ptr = nullptr;
};



template<typename T>
struct is_field_tile_set : std::false_type
{
};


template<typename GL, typename Arr, typename PQ>
struct is_field_tile_set<basic::FieldTileSet<GL, Arr, PQ>> : std::true_type
{
};

template<typename GL, typename Arr, typename PQ>
struct is_field_tile_set<FieldTileSet<GL, Arr, PQ>> : std::true_type
{
};

template<typename GL, typename Arr, typename PQ>
struct is_field_tile_set<GridTileSet<GL, Arr, PQ>> : std::true_type
{
};

template<typename T>
auto static constexpr is_field_tile_set_v = is_field_tile_set<T>::value;



template<typename GridLayout_t, typename Grid_t, typename Field_t>
inline std::ostream& operator<<(std::ostream& out,
                                FieldTileSet<GridLayout_t, Grid_t, Field_t> const& ts)
{
    for (auto const& tile : ts())
        out << tile();

    return out;
}


template<typename GridLayout_t, typename Grid_t, typename Field_t>
inline auto sum_field(FieldTileSet<GridLayout_t, Grid_t, Field_t> const& ts)
{
    return sum_from(ts(), [](auto const& tile) { return sum(tile()); });
}

template<typename GridLayout_t, typename Grid_t, typename Field_t>
inline auto sum_not_nan(FieldTileSet<GridLayout_t, Grid_t, Field_t> const& ts)
{
    return sum_from(ts(), [](auto const& tile) {
        typename Grid_t::value_type s = 0;
        for (auto const& v : tile())
            if (!std::isnan(v))
                s += v;
        return s;
    });
}


template<typename GridLayout_t, typename Grid_t, typename Field_t>
void check_field(basic::FieldTileSet<GridLayout_t, Grid_t, Field_t> const& f,
                 GridLayout_t const& layout)
{
#if PHARE_DEBUG
    auto const domainBox = layout.AMRBoxFor(f);

    for (auto& tile : f())
        if (auto const overlap = **tile * domainBox)
            for (auto const& bix : tile.layout().AMRToLocal(*overlap))
            {
                assert(not std::isnan(tile()(bix)));
            }
#endif // PHARE_DEBUG
}

template<typename GridLayout_t, typename Grid_t, typename Field_t>
void check_field(basic::FieldTileSet<GridLayout_t, Grid_t, Field_t> const& f)
{
#if PHARE_DEBUG
    for (auto& tile : f())
        check_field(tile());
#endif // PHARE_DEBUG
}



} // namespace PHARE::core



#endif // PHARE_CORE_DATA_GRID_GRID_TILES_HPP
