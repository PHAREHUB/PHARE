#ifndef PHARE_TILE_SET_HPP
#define PHARE_TILE_SET_HPP


#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "tile_set_mapper.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <tuple>


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

    NO_DISCARD auto& operator[](std::size_t const i) _PHARE_ALL_FN_ { return tiles_[i]; }
    NO_DISCARD auto& operator[](std::size_t const i) const _PHARE_ALL_FN_ { return tiles_[i]; }

    NO_DISCARD auto& operator()() { return tiles_; }
    NO_DISCARD auto& operator()() const { return tiles_; }



    template<typename... Index>
    NO_DISCARD auto at(Index... indexes) _PHARE_ALL_FN_
    {
        auto const& idx = cell_idx(indexes...);
        assert(idx < cells_.size());
        return cells_[idx];
    }
    template<typename... Index>
    NO_DISCARD auto at(Index... indexes) const _PHARE_ALL_FN_
    {
        auto const& idx = cell_idx(indexes...);
#if PHARE_DEBUG
        if (idx >= cells_.size())
        {
            PHARE_LOG_LINE_SS("out of bounds idx for size " << std::to_string(idx) << " : "
                                                            << std::to_string(cells_.size()));
            throw std::runtime_error("out of bounds idx for size " + std::to_string(idx) + " : "
                                     + std::to_string(cells_.size()));
        }
#endif

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
    auto cell_idx(auto const&... ijk) const _PHARE_ALL_FN_
    {
        return NdArrayViewer<dimension>::idx(cells_shape_, ijk...);
    }

    Box_t box_;
    Span_t tiles_;
    CellSpan_t cells_;
    std::array<std::uint32_t, dimension> cells_shape_;
};


template<std::size_t dim>
struct TileSetter
{
    using Box_t = Box<int, dim>;

    Box_t box;
    std::size_t ghosts;
};


template<typename Tile, auto alloc_mode = AllocatorMode::CPU>
class TileSet
{
    using This                    = TileSet<Tile, alloc_mode>;
    bool static constexpr c_order = true;

    template<typename T>
    using nd_array_t = NdArrayVector<Tile::dimension, T, c_order, alloc_mode>;

    using Point_t = Point<std::uint32_t, Tile::dimension>;

    TileSet() = default;

    template<typename, auto>
    friend class TileSet;

public:
    template<typename T>
    using vector_t = typename PHARE::Vector<T, alloc_mode>::vector_t;

    using value_type                = Tile;
    using Box_t                     = Box<int, Tile::dimension>;
    static auto constexpr dimension = Tile::dimension;


    TileSet(Box_t const& box, nd_array_t<Tile*> const& cells, vector_t<Tile> const& tiles)
        : box_{box}
        , cells_{cells}
        , tiles_{tiles}
    {
    }

    template<typename... Args>
    TileSet(TileSetter<Tile::dimension> const& setter, Args&&... args)
        : box_{setter.box}
        , cells_{grow(setter.box, setter.ghosts).shape().template toArray<std::uint32_t>()}
    {
        make_tiles_(args...);
        tag_cells_(setter);
        box_ = grow(setter.box, setter.ghosts); // ghost box now!
    }



    template<typename... Args>
    TileSet(Box_t const& box, Args&&... args)
        : box_{box}
        , cells_{box.shape().template toArray<std::uint32_t>()}
    {
        make_tiles_(args...);
        tag_cells_();
    }


    TileSet(TileSet const& that)
        : box_{that.box_}
        , cells_{that.cells_}
        , tiles_{that.tiles_}
    {
        tag_cells_();
    }


    TileSet copy(TileSetter<Tile::dimension> const& setter) const
    {
        TileSet ret{box_, cells_, tiles_};
        ret.tag_cells_(setter);
        return ret;
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

    auto static make_from(auto const& accessor, TileSetter<Tile::dimension> const& setter,
                          auto& that, auto&&... args)
    {
        This ts{setter.box,
                nd_array_t<Tile*>{
                    grow(setter.box, setter.ghosts).shape().template toArray<std::uint32_t>()},
                {}};



        // ts.box_   = that.box_;
        // ts.cells_ = {ts.cells_.shape()};
        for (auto& tile : that())
            ts.tiles_.emplace_back(accessor(tile), args...);
        ts.tag_cells_(setter);
        return ts;
    }

    // auto static make_from(auto const& that, auto&&... args)
    // {
    //     This ts{};
    //     ts.box_   = that.box_;
    //     ts.cells_ = {ts.cells_.shape()};
    //     for (auto const& tile : that())
    //         ts.tiles_.emplace_back(tile, args...);
    //     ts.tag_cells_();
    //     return ts;
    // }


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

        if constexpr (dimension == 3)
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


    template<typename... Args>
    void make_tiles_(Args&&... args)
    {
        tile_set_make_tiles(*this, args...);
    }

    void tag_cells_(auto const& setter)
    {
        // use setter.box (domain box) explicitly — box_ may already be the ghost box
        // when called from copy(), so we cannot rely on it here
        auto const& domain_box = setter.box;
        auto const ghost_box   = grow(domain_box, setter.ghosts);

        for (auto& tile : tiles_)
            for (auto const& cell : tile)
                cells_(cell - ghost_box.lower) = &tile;

        for (auto& tile : tiles_)
            if (auto const ghost_tile = grow(tile, setter.ghosts);
                ghost_tile * domain_box != ghost_tile)
                for (auto const& cell : ghost_tile)
                    if (!isIn(cell, domain_box))
                    {
                        auto clamped = cell;
                        for (std::size_t i = 0; i < dimension; ++i)
                            clamped[i] = std::min(std::max(cell[i], domain_box.lower[i]),
                                                  domain_box.upper[i]);
                        cells_(cell - ghost_box.lower) = cells_(clamped - ghost_box.lower);
                    }
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
    vector_t<Tile> tiles_{};
};


template<typename F, typename Tile, auto alloc_mode = AllocatorMode::CPU>
auto& update_from(F f, TileSet<Tile, alloc_mode>& in)
{
    for (std::size_t i = 0; i < in.size(); ++i)
        in.data()[i] = f(i);
    return in;
}

// template<auto alloc_mode = AllocatorMode::CPU, typename F, typename Tile, auto am1,
//          typename... Args>
// auto generate_from(F f, TileSet<Tile, am1>& in, Args&&... args)
// {
//     using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
//     TileSet<value_type, alloc_mode> ret{in.box(), args...};
//     assert(in.size() == ret.size());

//     for (std::size_t i = 0; i < in.size(); ++i)
//         ret.data()[i] = f(i);
//     return ret;
// }

// template<auto alloc_mode = AllocatorMode::CPU, typename F, typename Tile, auto am1,
//          typename... Args>
// auto generate_from(F f, TileSet<Tile, am1> const& in, Args&&... args)
// {
//     using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
//     TileSet<value_type, alloc_mode> ret{in.box(), args...};
//     assert(in.size() == ret.size());

//     for (std::size_t i = 0; i < in.size(); ++i)
//         ret.data()[i] = f(i);
//     return ret;
// }

template<auto alloc_mode = AllocatorMode::CPU, std::size_t dim, typename F, typename Tile, auto am1,
         typename... Args>
auto generate_from(TileSetter<dim> const& setter, F f, TileSet<Tile, am1> const& in, Args&&... args)
{
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    TileSet<value_type, alloc_mode> ret{setter, args...};
    assert(in.size() == ret.size());

    for (std::size_t i = 0; i < in.size(); ++i)
        ret.data()[i] = f(i);
    return ret;
}

} // namespace PHARE::core


#endif
