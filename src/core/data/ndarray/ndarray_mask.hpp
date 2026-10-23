#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_MASK_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_MASK_HPP


#include <array>
#include <tuple>
#include <vector>
#include <cstdint>
#include <numeric>
#include <iostream>
#include <stdexcept>


#include "core/def.hpp"
#include "core/vector.hpp"
#include "core/utilities/types.hpp"

#include "core/data/ndarray/ndarray_base.hpp"

namespace PHARE::core
{

template<typename Array, typename Mask>
class MaskedView
{
public:
    static auto constexpr dimension = Array::dimension;
    using DataType                  = typename Array::type;
    using data_type                 = typename Array::type;

    MaskedView(Array& array, Mask const& mask)
        : array_{array}
        , shape_{array.shape()}
        , mask_{mask}
    {
    }

    MaskedView(Array& array, Mask&& mask)
        : array_{array}
        , shape_{array.shape()}
        , mask_{std::move(mask)}
    {
    }

    template<typename... Indexes>
    NO_DISCARD DataType const& operator()(Indexes... indexes) const
    {
        return NdArrayViewer<dimension, true>::at(array_.data(), shape_, indexes...);
    }

    template<typename... Indexes>
    NO_DISCARD DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<MaskedView const&>(*this)(indexes...));
    }

    NO_DISCARD auto operator=(data_type value) { mask_.fill(array_, value); }

    NO_DISCARD auto xstart() const { return mask_.min(); }
    NO_DISCARD auto xend() const { return shape_[0] - 1 - mask_.max(); }

    NO_DISCARD auto ystart() const { return mask_.min(); }
    NO_DISCARD auto yend() const { return shape_[1] - 1 - mask_.max(); }

    NO_DISCARD auto zstart() const { return mask_.min(); }
    NO_DISCARD auto zend() const { return shape_[2] - 1 - mask_.max(); }


private:
    Array& array_;
    std::array<std::uint32_t, dimension> shape_;
    Mask const& mask_;
};




class NdArrayMask
{
public:
    NdArrayMask(std::size_t min, std::size_t max)
        : min_{min}
        , max_{max}
    {
    }

    NdArrayMask(std::size_t width)
        : min_{width}
        , max_{width}
    {
    }

    template<typename Array>
    void fill(Array& array, typename Array::type val) const
    {
        if constexpr (Array::dimension == 1)
            fill1D(array, val);

        else if constexpr (Array::dimension == 2)
            fill2D(array, val);

        else if constexpr (Array::dimension == 3)
            fill3D(array, val);
    }

    template<typename Array>
    void fill1D(Array& array, typename Array::type val) const
    {
        auto shape = array.shape();

        for (std::size_t i = min_; i <= max_; ++i)
            array(i) = val;

        for (std::size_t i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            array(i) = val;
    }

    template<typename Array>
    void fill2D(Array& array, typename Array::type val) const
    {
        auto shape = array.shape();

        // left border
        for (std::size_t i = min_; i <= max_; ++i)
            for (std::size_t j = min_; j <= shape[1] - 1 - max_; ++j)
                array(i, j) = val;

        // right border
        for (std::size_t i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            for (std::size_t j = min_; j <= shape[1] - 1 - max_; ++j)
                array(i, j) = val;


        for (std::size_t i = min_; i <= shape[0] - 1 - min_; ++i)
        {
            // bottom border
            for (std::size_t j = min_; j <= max_; ++j)
                array(i, j) = val;

            // top border
            for (std::size_t j = shape[1] - 1 - max_; j <= shape[1] - 1 - min_; ++j)
                array(i, j) = val;
        }
    }


    template<typename Array>
    void fill3D(Array& array, typename Array::type val) const
    {
        auto shape = array.shape();

        // left border
        for (auto i = min_; i <= shape[0] - 1 - max_; ++i)
            for (auto j = min_; j <= shape[1] - 1 - max_; ++j)
                for (auto k = min_; k <= max_; ++k)
                    array(i, j, k) = val;

        // // right border
        for (auto i = min_; i <= shape[0] - 1 - max_; ++i)
            for (auto j = min_; j <= shape[1] - 1 - max_; ++j)
                for (auto k = shape[2] - 1 - max_; k <= shape[2] - 1 - min_; ++k)
                    array(i, j, k) = val;

        for (auto i = min_; i <= shape[0] - 1 - min_; ++i)
        {
            // bottom border
            for (auto j = min_; j <= max_; ++j)
                for (auto k = min_; k <= shape[2] - 1 - min_; ++k)
                    array(i, j, k) = val;

            // top border
            for (auto j = shape[1] - 1 - max_; j <= shape[1] - 1 - min_; ++j)
                for (auto k = min_; k <= shape[2] - 1 - min_; ++k)
                    array(i, j, k) = val;
        }

        // front
        for (auto i = min_; i <= max_; ++i)
            for (auto j = min_; j <= shape[1] - 1 - max_; ++j)
                for (auto k = min_; k <= shape[2] - 1 - min_; ++k)
                    array(i, j, k) = val;

        // back
        for (auto i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            for (auto j = min_; j <= shape[1] - 1 - max_; ++j)
                for (auto k = min_; k <= shape[2] - 1 - min_; ++k)
                    array(i, j, k) = val;
    }

    template<typename Array>
    NO_DISCARD auto nCells(Array const& array)
    {
        auto shape = array.shape();

        std::size_t cells = 0;

        for (auto i = min_; i <= max_; ++i)
        {
            if constexpr (Array::dimension == 1)
                cells += 2;

            if constexpr (Array::dimension == 2)
                cells += (shape[0] - (i * 2) - 2) * 2 + (shape[1] - (i * 2) - 2) * 2 + 4;

            if constexpr (Array::dimension == 3)
            {
                auto [x, y, z] = shape;
                x -= i * 2;
                y -= i * 2;
                z -= i * 2;
                cells += (x * y * 2) + (y * (z - 2) * 2) + ((z - 2) * (x - 2) * 2);
            }
        }
        return cells;
    }


    NO_DISCARD auto min() const { return min_; };
    NO_DISCARD auto max() const { return max_; };

private:
    std::size_t min_ = 0, max_ = 0;
};




template<typename Array, typename Mask>
void operator>>(MaskedView<Array, Mask>&& inner, MaskedView<Array, Mask>&& outer)
{
    using MaskedView_t = MaskedView<Array, Mask>;

    assert(inner.xstart() > outer.xstart() and inner.xend() < outer.xend());

    if constexpr (MaskedView_t::dimension == 1)
    {
        outer(outer.xstart()) = inner(inner.xstart());
        outer(outer.xend())   = inner(inner.xend());
    }

    if constexpr (MaskedView_t::dimension > 1)
    {
        assert(inner.ystart() > outer.ystart() and inner.yend() < outer.yend());
    }

    if constexpr (MaskedView_t::dimension == 2)
    {
        for (auto ix = inner.xstart(); ix <= inner.xend(); ++ix)
        {
            outer(ix, outer.ystart()) = inner(ix, inner.ystart()); // bottom
            outer(ix, outer.yend())   = inner(ix, inner.yend());   // top
        }

        for (auto iy = inner.ystart(); iy <= inner.yend(); ++iy)
        {
            outer(outer.xstart(), iy) = inner(inner.xstart(), iy); // left
            outer(outer.xend(), iy)   = inner(inner.xend(), iy);   // right
        }

        // bottom left
        for (auto ix = outer.xstart(); ix < inner.xstart(); ++ix)
            outer(ix, outer.ystart()) = inner(inner.xstart(), inner.ystart());

        for (auto iy = outer.ystart(); iy < inner.ystart(); ++iy)
            outer(outer.xstart(), iy) = inner(inner.xstart(), inner.ystart());


        // top left
        for (auto ix = outer.xstart(); ix < inner.xstart(); ++ix)
            outer(ix, outer.yend()) = inner(inner.xstart(), inner.yend());

        for (auto iy = outer.yend(); iy > inner.yend(); --iy)
            outer(outer.xstart(), iy) = inner(inner.xstart(), inner.yend());

        // top right
        for (auto ix = outer.xend(); ix > inner.xend(); --ix)
            outer(ix, outer.yend()) = inner(inner.xend(), inner.yend());

        for (auto iy = outer.yend(); iy > inner.yend(); --iy)
            outer(outer.xend(), iy) = inner(inner.xend(), inner.yend());


        // bottom right
        for (auto ix = outer.xend(); ix > inner.xend(); --ix)
            outer(ix, outer.ystart()) = inner(inner.xend(), inner.ystart());

        for (auto iy = outer.ystart(); iy < inner.ystart(); ++iy)
            outer(outer.xend(), iy) = inner(inner.xend(), inner.ystart());
    }


    if constexpr (MaskedView_t::dimension == 3)
    {
        assert(inner.zstart() > outer.zstart() and inner.zend() < outer.zend());

        for (auto i = inner.xstart(); i <= inner.xend(); ++i)
        {
            for (auto k = inner.zstart(); k <= inner.zend(); ++k)
            {
                outer(i, outer.ystart(), k) = inner(i, inner.ystart(), k);
                outer(i, outer.yend(), k)   = inner(i, inner.yend(), k);
            }
            for (auto j = inner.ystart(); j <= inner.yend(); ++j)
            {
                outer(i, j, outer.zstart()) = inner(i, j, inner.zstart());
                outer(i, j, outer.zend())   = inner(i, j, inner.zend());
            }
        }

        for (auto j = inner.ystart(); j <= inner.yend(); ++j)
        {
            for (auto k = inner.zstart(); k <= inner.zend(); ++k)
            {
                outer(outer.xstart(), j, k) = inner(inner.xstart(), j, k);
                outer(outer.xend(), j, k)   = inner(inner.xend(), j, k);
            }
        }

        for (auto k = inner.zstart(); k <= inner.zend(); ++k)
        {
            outer(outer.xstart(), outer.ystart(), k) = inner(inner.xstart(), inner.ystart(), k);
            outer(outer.xstart(), outer.yend(), k)   = inner(inner.xstart(), inner.yend(), k);

            outer(outer.xend(), outer.ystart(), k) = inner(inner.xend(), inner.ystart(), k);
            outer(outer.xend(), outer.yend(), k)   = inner(inner.xend(), inner.yend(), k);
        }

        for (auto j = inner.ystart(); j <= inner.yend(); ++j)
        {
            outer(outer.xstart(), j, outer.zstart()) = inner(inner.xstart(), j, inner.zstart());
            outer(outer.xstart(), j, outer.zend())   = inner(inner.xstart(), j, inner.zend());

            outer(outer.xend(), j, outer.zstart()) = inner(inner.xend(), j, inner.zstart());
            outer(outer.xend(), j, outer.zend())   = inner(inner.xend(), j, inner.zend());
        }

        for (auto i = inner.xstart(); i <= inner.xend(); ++i)
        {
            outer(i, outer.ystart(), outer.zstart()) = inner(i, inner.ystart(), inner.zstart());
            outer(i, outer.ystart(), outer.zend())   = inner(i, inner.ystart(), inner.zend());

            outer(i, outer.yend(), outer.zstart()) = inner(i, inner.yend(), inner.zstart());
            outer(i, outer.yend(), outer.zend())   = inner(i, inner.yend(), inner.zend());
        }

        auto corner = [&](auto xouter, auto xinner) {
            outer(xouter, outer.ystart(), outer.zstart())
                = inner(xinner, inner.ystart(), inner.zstart());

            outer(xouter, outer.ystart(), outer.zend())
                = inner(xinner, inner.ystart(), inner.zend());

            outer(xouter, outer.yend(), outer.zstart())
                = inner(xinner, inner.yend(), inner.zstart());
            outer(xouter, outer.yend(), outer.zend()) = inner(xinner, inner.yend(), inner.zend());
        };

        corner(outer.xstart(), inner.xstart());
        corner(outer.xend(), inner.xend());
    }
}

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_MASK_HPP
