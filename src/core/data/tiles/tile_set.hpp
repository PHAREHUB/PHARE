#ifndef PHARE_TILE_SET_HPP
#define PHARE_TILE_SET_HPP

#include "core/utilities/box/box.hpp"
#include "core/utilities/types.hpp"

#include <iostream>
#include <array>
#include <tuple>
#include <string>

namespace PHARE::core
{
template<typename Tile>
class TileSet
{
public:
    static auto constexpr dimension = Tile::dimension;

    TileSet(Box<int, dimension> const& box, std::array<int, dimension> const& tile_size)
        : box_{box}
        , tile_size_{tile_size}
        , shape_{[&]() {
            std::array<int, dimension> s;
            auto bs = box.shape();
            for (auto i = 0u; i < dimension; ++i)
            {
                auto const div = (bs[i] + tile_size_[i] - 1) / tile_size_[i];
                s[i]           = div;
            }
            return s;
        }()}
        , tiles_(product(shape_))
    {
        for (auto idim = 0u; idim < dimension; ++idim)
        {
            if (box_.shape()[idim] < tile_size_[idim])
            {
                throw std::runtime_error("tile size larger than box size in dimension "
                                         + std::to_string(idim));
            }
        }

        auto const size_me = [&](auto dim, auto idx) {
            if (idx == shape_[dim] - 1)
            {
                auto const remain = box_.shape()[dim] % tile_size_[dim];
                return (remain == 0) ? tile_size_[dim] : remain;
            }
            else
                return tile_size_[dim];
        };

        for (auto ix = 0u; ix < shape_[0]; ++ix)
        {
            if constexpr (dimension == 1)
            {
                // -1 because upper is included
                tiles_[ix].lower[0] = box.lower[0] + ix * tile_size_[0];
                tiles_[ix].upper[0] = tiles_[ix].lower[0] + size_me(0, ix) - 1;
            }
            else
            {
                for (auto iy = 0u; iy < shape_[1]; ++iy)
                {
                    if constexpr (dimension == 2)
                    {
                        auto const i       = ix * shape_[1] + iy;
                        tiles_[i].lower[0] = box.lower[0] + ix * tile_size_[0];
                        tiles_[i].upper[0] = tiles_[i].lower[0] + size_me(0, ix) - 1;
                        tiles_[i].lower[1] = box.lower[1] + iy * tile_size_[1];
                        tiles_[i].upper[1] = tiles_[i].lower[1] + size_me(1, iy) - 1;
                    }
                    else
                    {
                        for (auto iz = 0u; iz < shape_[2]; ++iz)
                        {
                            auto const i       = ix * shape_[1] * shape_[2] + shape_[2] * iy + iz;
                            tiles_[i].lower[0] = box.lower[0] + ix * tile_size_[0];
                            tiles_[i].upper[0] = tiles_[i].lower[0] + size_me(0, ix) - 1;
                            tiles_[i].lower[1] = box.lower[1] + iy * tile_size_[1];
                            tiles_[i].upper[1] = tiles_[i].lower[1] + size_me(1, iy) - 1;
                            tiles_[i].lower[2] = box.lower[2] + iz * tile_size_[2];
                            tiles_[i].upper[2] = tiles_[i].lower[2] + size_me(2, iz) - 1;
                        }
                    }
                }
            }
        }
    }


    auto export_overlaped_with(Box<int, dimension> const& box) const
    {
        std::vector<std::pair<bool, Tile>> overlaped;
        for (auto const& tile : tiles_)
        {
            auto overlap = box * Box<int, dimension>{tile.lower, tile.upper};
            if (overlap)
            {
                auto complete_overlap = (*overlap).size() == tile.size();
                overlaped.emplace_back(complete_overlap, tile);
            }
        }
        return overlaped;
    }

    auto shape() const { return shape_; }
    auto size() const { return tiles_.size(); }

    auto begin() { return tiles_.begin(); }
    auto begin() const { return tiles_.begin(); }

    auto end() { return tiles_.end(); }
    auto end() const { return tiles_.end(); }

    auto& operator[](std::size_t i) { return tiles_[i]; }
    auto const& operator[](std::size_t i) const { return tiles_[i]; }

private:
    Box<int, dimension> box_;
    std::array<int, dimension> tile_size_;
    std::array<int, dimension> shape_;
    std::vector<Tile> tiles_;
};
} // namespace PHARE::core


#endif
