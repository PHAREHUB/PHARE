#ifndef PHARE_TILE_SET_HPP
#define PHARE_TILE_SET_HPP

#include "core/def.hpp"
#include "core/data/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "tile_set_mapper.hpp"

#include <array>
#include <tuple>
#include <vector>

namespace PHARE::core
{

template<typename Tile, typename Span_t = Span<Tile>, typename CellSpan_t = Span<Tile*>>
class TileSetView
{
public:
    using value_type                = Tile;
    using Box_t                     = Box<int, Tile::dimension>;
    static auto constexpr dimension = Tile::dimension;

    TileSetView(Box_t const& box, Tile* tiles, std::size_t tile_nbr, Tile** cells,
                std::array<std::uint32_t, dimension> const& nbr_cells)
        : box_{box}
        , tiles_{tiles, tile_nbr}
        , cells_{cells, product(nbr_cells)}
        , cells_shape_{nbr_cells}
    {
    }

    TileSetView(TileSetView const&)            = default;
    TileSetView(TileSetView&&)                 = default;
    TileSetView& operator=(TileSetView const&) = default;
    TileSetView& operator=(TileSetView&&)      = default;


    NO_DISCARD auto box() const { return box_; }
    NO_DISCARD auto size() const { return tiles_.size(); }

    NO_DISCARD auto begin() { return tiles_.begin(); }
    NO_DISCARD auto begin() const { return tiles_.begin(); }

    NO_DISCARD auto end() { return tiles_.end(); }
    NO_DISCARD auto end() const { return tiles_.end(); }

    NO_DISCARD auto data() { return tiles_.data(); }
    NO_DISCARD auto data() const { return tiles_.data(); }

    NO_DISCARD auto& operator[](std::size_t const i) { return tiles_[i]; }
    NO_DISCARD auto& operator[](std::size_t const i) const { return tiles_[i]; }

    NO_DISCARD auto& operator()() { return tiles_; }
    NO_DISCARD auto& operator()() const { return tiles_; }

    template<typename... Index>
    NO_DISCARD auto at(Index... indexes)
    {
        auto const& idx = cell_idx(indexes...);
        assert(idx < cells_.size());
        return cells_[idx];
    }
    template<typename... Index>
    NO_DISCARD auto at(Index... indexes) const
    {
        auto const& idx = cell_idx(indexes...);
        return cells_[idx];
    }

    template<typename Tile_t>
    void reset(std::size_t i, Tile_t& tile)
    {
        tiles_[i].reset(tile);
    }

    template<typename TileSet_t>
    void reset(TileSet_t& tileset)
    {
        for (std::size_t i = 0; i < tileset.size(); ++i)
            reset(i, tileset[i]);
    }

    void reset() { tiles_.ptr = nullptr; }

protected:
    auto cell_idx(auto const&... ijk) const
    {
        return NdArrayViewer<dimension>::idx(cells_shape_, ijk...);
    }

    Box_t box_;
    Span_t tiles_;
    CellSpan_t cells_;
    std::array<std::uint32_t, dimension> cells_shape_;
};

template<typename Tile, auto alloc_mode = AllocatorMode::CPU>
class TileSet
{
    bool static constexpr c_order = true;
    using Point_t                 = Point<std::uint32_t, Tile::dimension>;

    TileSet() = default;

    template<typename, auto>
    friend class TileSet;

protected:
    template<typename T>
    using nd_array_t = NdArrayVector<Tile::dimension, T, c_order, alloc_mode>;

public:
    using value_type                = Tile;
    using Box_t                     = Box<int, Tile::dimension>;
    static auto constexpr dimension = Tile::dimension;

    TileSet(Box_t const& box, nd_array_t<Tile*> const& cells, std::vector<Tile> const& tiles)
        : box_{box}
        , cells_{cells}
        , tiles_{tiles}
    {
    }

    template<typename... Args>
    TileSet(Box_t const& box, Args&&... args)
        : box_{box}
        , cells_{box.shape().template toArray<std::uint32_t>()}
    {
        make_tiles_(TilingOptions::from_env(), args...);
        tag_cells_();
    }

    TileSet(TileSet const& that)
        : box_{that.box_}
        , cells_{that.cells_}
        , tiles_{that.tiles_}
    {
        tag_cells_();
    }

    TileSet(TileSet&&) = default;

    TileSet& operator=(TileSet const& that)
    {
        box_   = that.box_;
        cells_ = that.cells_;
        tiles_ = that.tiles_;
        tag_cells_();
        return *this;
    }

    NO_DISCARD auto box() const { return box_; }
    NO_DISCARD auto size() const { return tiles_.size(); }

    NO_DISCARD auto begin() { return tiles_.begin(); }
    NO_DISCARD auto begin() const { return tiles_.begin(); }

    NO_DISCARD auto end() { return tiles_.end(); }
    NO_DISCARD auto end() const { return tiles_.end(); }

    NO_DISCARD auto data() { return tiles_.data(); }
    NO_DISCARD auto data() const { return tiles_.data(); }

    NO_DISCARD auto& operator[](std::size_t i) { return tiles_[i]; }
    NO_DISCARD auto const& operator[](std::size_t i) const { return tiles_[i]; }

    NO_DISCARD auto& operator()() { return tiles_; }
    NO_DISCARD auto& operator()() const { return tiles_; }

    NO_DISCARD auto& at(Point<int, dimension> const& amr_point)
    {
        return cells_((amr_point - box_.lower).as_unsigned());
    }
    NO_DISCARD auto& at(Point<int, dimension> const& amr_point) const
    {
        return cells_((amr_point - box_.lower).as_unsigned());
    }

    template<typename... Index>
    NO_DISCARD auto& at(Index... indexes)
    {
        return cells_(indexes...);
    }
    template<typename... Index>
    NO_DISCARD auto& at(Index... indexes) const
    {
        return cells_(indexes...);
    }

    template<typename View_t = Tile>
    auto make_view()
    {
        return TileSetView<View_t>{box_, tiles_.data(), tiles_.size(), cells_.data(),
                                   cells_.shape()};
    }
    template<typename View_t = Tile>
    auto make_view() const
    {
        return TileSetView<View_t>{box_, const_cast<Tile*>(tiles_.data()), tiles_.size(),
                                   const_cast<Tile**>(cells_.data()), cells_.shape()};
    }

    auto as(auto&& a, auto&&... args)
    {
        return a(box_, tiles_.data(), tiles_.size(), cells_.data(), cells_.shape(), args...);
    }

    void static build_links(TileSet& tile_set)
    {
        auto const box_shape = tile_set.box_.shape().as_unsigned();

        if constexpr (dimension == 1)
        {
            for (auto xi = 0u; xi < box_shape[0];)
            {
                Point const xip{xi};
                auto const xitile = tile_set.cells_(xip);
                link_dim1(tile_set, box_shape, xip);
                xi += xitile->shape().as_unsigned()[0];
            }
        }
        else if constexpr (dimension == 2)
        {
            for (auto xi = 0u; xi < box_shape[0];)
            {
                Point const xip{xi, 0u};
                auto const xitile       = tile_set.cells_(xip);
                auto const xitile_shape = xitile->shape().as_unsigned();

                for (auto yi = 0u; yi < box_shape[1];)
                {
                    Point const yip{xi, yi};
                    auto const yitile = tile_set.cells_(yip);
                    link_dim2(tile_set, box_shape, yip);
                    yi += yitile->shape().as_unsigned()[1];
                }
                xi += xitile_shape[0];
            }
        }
        else if constexpr (dimension == 3)
        {
            // links[0] = {0, 0, 1};
            // links[1] = {0, 1, 0};
            // links[2] = {0, 1, 1};
            // links[3] = {1, 0, 0};
            // links[4] = {1, 0, 1};
            // links[5] = {1, 1, 0};
            // links[6] = {1, 1, 1};

            for (auto xi = 0u; xi < box_shape[0];)
            {
                Point const xip{xi, 0u, 0u};
                auto const xitile       = tile_set.cells_(xip);
                auto const xitile_shape = xitile->shape().as_unsigned();

                for (auto yi = 0u; yi < box_shape[1];)
                {
                    Point const yip{xi, yi, 0u};
                    auto const yitile       = tile_set.cells_(yip);
                    auto const yitile_shape = yitile->shape().as_unsigned();

                    for (auto zi = 0u; zi < box_shape[2];)
                    {
                        Point const zip{xi, yi, zi};
                        auto const zitile       = tile_set.cells_(zip);
                        auto const zitile_shape = zitile->shape().as_unsigned();
                        link_dim3(tile_set, box_shape, zip);
                        zi += zitile_shape[2];
                    }
                    yi += yitile_shape[1];
                }
                xi += xitile_shape[0];
            }
        }
    }

private:
    template<typename... Args>
    void static _link(Args&&... args)
    {
        auto const& [tile_set, tile, box_shape, point, idx] = std::forward_as_tuple(args...);
        if (for_N_all<dimension>([&](auto i) { return point[i] < box_shape[i]; }))
            tile->link(idx) = tile_set.cells_(point);
    }

    template<typename... Args>
    void static link_dim1(Args&&... args)
    {
        auto const& [tile_set, box_shape, point] = std::forward_as_tuple(args...);
        auto tile                                = tile_set.cells_(point);
        auto const tile_shape                    = tile->shape().as_unsigned();

        _link(tile_set, tile, box_shape, point + tile_shape, 0);
    }

    template<typename... Args>
    void static link_dim2(Args&&... args)
    {
        auto const& [tile_set, box_shape, point] = std::forward_as_tuple(args...);
        auto tile                                = tile_set.cells_(point);
        auto const tile_shape                    = tile->shape().as_unsigned();

        { // +Y (fastest/last axis, mirrors link_dim3's slot 0)
            auto link0 = point;
            link0[1] += tile_shape[1];
            _link(tile_set, tile, box_shape, link0, 0);
        }
        { // +X (mirrors link_dim3's slot 1/3 role as the outer axis)
            auto link1 = point;
            link1[0] += tile_shape[0];
            _link(tile_set, tile, box_shape, link1, 1);
        }
        // +X+Y corner
        _link(tile_set, tile, box_shape, point + tile_shape, 2);
    }

    template<typename... Args>
    void static link_dim3(Args&&... args)
    {
        auto const& [tile_set, box_shape, point] = std::forward_as_tuple(args...);
        auto tile                                = tile_set.cells_(point);
        auto const tile_shape                    = tile->shape().as_unsigned();

        {
            auto link0 = point;
            link0[2] += tile_shape[2];
            _link(tile_set, tile, box_shape, link0, 0);
        }
        {
            auto link1 = point;
            link1[1] += tile_shape[1];
            _link(tile_set, tile, box_shape, link1, 1);
        }
        {
            auto link2 = point;
            link2[1] += tile_shape[1];
            link2[2] += tile_shape[2];
            _link(tile_set, tile, box_shape, link2, 2);
        }
        {
            auto link3 = point;
            link3[0] += tile_shape[0];
            _link(tile_set, tile, box_shape, link3, 3);
        }
        {
            auto link4 = point;
            link4[0] += tile_shape[0];
            link4[2] += tile_shape[2];
            _link(tile_set, tile, box_shape, link4, 4);
        }
        {
            auto link5 = point;
            link5[0] += tile_shape[0];
            link5[1] += tile_shape[1];
            _link(tile_set, tile, box_shape, link5, 5);
        }

        _link(tile_set, tile, box_shape, point + tile_shape, 6);
    }

protected:
    template<typename... Args>
    void make_tiles_(TilingOptions const& tiling, Args&&... args)
    {
        tile_set_make_tiles(*this, tiling, args...);
    }

    //! store the pointer to the tile associated with each cell
    void tag_cells_()
    {
        for (auto& tile : tiles_)         // need to substract box lower to get
            for (auto const& cell : tile) // the local index of that cell in the NdArray
                cells_(cell - box_.lower) = &tile;
    }

    Box_t box_;
    nd_array_t<Tile*> cells_;
    std::vector<Tile> tiles_{};
};


// tiles are built from the domain box, but `.at()` covers the ghost box: each ghost cell maps to
// the tile owning the clamped (nearest) domain cell, so every ghost cell has a single owner tile
template<typename Tile, auto alloc_mode = AllocatorMode::CPU>
class GhostClampedTileSet : public TileSet<Tile, alloc_mode>
{
    using Super = TileSet<Tile, alloc_mode>;

    template<typename T>
    using nd_array_t = Super::template nd_array_t<T>;

    template<typename, auto>
    friend class GhostClampedTileSet;

    struct untiled_t
    {
    };

    // `{}` tiles: keeps TileSet's variadic tiling constructor out of overload resolution
    GhostClampedTileSet(untiled_t, Super::Box_t const& domain_box, std::size_t const ghosts)
        : Super{domain_box, ghost_cells_(domain_box, ghosts), {}}
        , domain_box_{domain_box}
        , ghosts_{ghosts}
    {
    }

public:
    using Box_t = Super::Box_t;
    using Super::dimension;

    template<typename... Args>
    GhostClampedTileSet(Box_t const& domain_box, std::size_t const ghosts, Args&&... args)
        : GhostClampedTileSet{untiled_t{}, domain_box, ghosts}
    {
        this->make_tiles_(TilingOptions::from_env(), args...); // tiles the domain box
        tag_cells_();
    }

    GhostClampedTileSet(GhostClampedTileSet const& that)
        : Super{that.box_, that.cells_, that.tiles_}
        , domain_box_{that.domain_box_}
        , ghosts_{that.ghosts_}
    {
        tag_cells_();
    }

    GhostClampedTileSet(GhostClampedTileSet&&) = default;

    GhostClampedTileSet& operator=(GhostClampedTileSet const& that)
    {
        this->box_   = that.box_;
        this->cells_ = that.cells_;
        this->tiles_ = that.tiles_;
        domain_box_  = that.domain_box_;
        ghosts_      = that.ghosts_;
        tag_cells_();
        return *this;
    }

    // one tile per tile of `that`, constructed from accessor(tile), args...
    auto static make_from(auto const& accessor, auto& that, auto&&... args)
    {
        GhostClampedTileSet ts{untiled_t{}, that.domain_box_, that.ghosts_};
        for (auto& tile : that())
            ts.tiles_.emplace_back(accessor(tile), args...);
        ts.tag_cells_();
        return ts;
    }

    NO_DISCARD auto& domain_box() const { return domain_box_; }
    NO_DISCARD auto ghosts() const { return ghosts_; }

private:
    auto static ghost_cells_(Box_t const& domain_box, std::size_t const ghosts)
    {
        return nd_array_t<Tile*>{
            grow(domain_box, ghosts).shape().template toArray<std::uint32_t>()};
    }

    void tag_cells_()
    {
        auto const ghost_box = grow(domain_box_, ghosts_);
        auto& cells          = this->cells_;

        for (auto& tile : this->tiles_)
            for (auto const& cell : tile)
                cells(cell - ghost_box.lower) = &tile;

        for (auto& tile : this->tiles_)
            if (auto const ghost_tile = grow(tile, ghosts_); ghost_tile * domain_box_ != ghost_tile)
                for (auto const& cell : ghost_tile)
                    if (!isIn(cell, domain_box_))
                    {
                        auto clamped = cell;
                        for (std::size_t i = 0; i < dimension; ++i)
                            clamped[i] = std::min(std::max(cell[i], domain_box_.lower[i]),
                                                  domain_box_.upper[i]);
                        cells(cell - ghost_box.lower) = cells(clamped - ghost_box.lower);
                    }

        this->box_ = ghost_box;
    }

    Box_t domain_box_;
    std::size_t ghosts_;
};


template<typename F, typename Tile, auto alloc_mode = AllocatorMode::CPU>
auto& update_from(F f, TileSet<Tile, alloc_mode>& in)
{
    for (std::size_t i = 0; i < in.size(); ++i)
        in.data()[i] = f(i);
    return in;
}

} // namespace PHARE::core


#endif
