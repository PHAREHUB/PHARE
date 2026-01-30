#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_HPP

#include "core/def.hpp"
#include <stdexcept>
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>
#include <tuple>
#include <numeric>
#include <iostream>


#include "core/utilities/types.hpp"


namespace PHARE::core
{
template<std::size_t dim, bool c_ordering = true>
struct NdArrayViewer
{
    using Id  = std::uint32_t;
    using Idx = Id const;

    template<typename NCells, template<typename, std::size_t> typename Indexes, typename Index>
    static inline std::uint32_t idx(NCells const& nCells, Indexes<Index, dim> const& indexes)
    {
        if constexpr (dim == 1)
            return idx(nCells, indexes[0]);

        else if constexpr (dim == 2)
            return idx(nCells, indexes[0], indexes[1]);

        else if constexpr (dim == 3)
            return idx(nCells, indexes[0], indexes[1], indexes[2]);
    }

    static inline std::uint32_t idx(auto const /*nCells*/, Idx i) { return i; }

    static inline std::uint32_t idx(auto const nCells, Idx i, Idx j)
    {
        if constexpr (c_ordering)
            return j + i * nCells[1];
        else
            return i + j * nCells[0];
    }
    static inline std::uint32_t idx(auto const nCells, Idx i, Idx j, Idx k)
    {
        if constexpr (c_ordering)
            return k + j * nCells[2] + i * nCells[1] * nCells[2];
        else
            return i + j * nCells[0] + k * nCells[1] * nCells[0];
    }

    template<template<typename, std::size_t> typename Indexes, typename Index>
    NO_DISCARD static inline auto& at(auto* data, auto const& nCells,
                                      Indexes<Index, dim> const& indexes)

    {
        auto const& i = idx(nCells, indexes);
        assert(i < product(nCells, std::uint32_t{1}));
        return data[i];
    }

    static inline auto& at(auto* data, auto const nCells, auto const... indexes)
    {
        auto const& i = idx(nCells, indexes...);
        assert(i < product(nCells, std::uint32_t{1}));
        return data[i];
    }
};




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



template<std::size_t dim, typename DataType = double, bool c_ordering = true>
class NdArrayView
{
public:
    static constexpr bool is_contiguous = 1;
    static std::size_t const dimension  = dim;
    using type                          = DataType;
    using pointer_type                  = DataType*;
    using viewer                        = NdArrayViewer<dim, c_ordering>;

    explicit NdArrayView(pointer_type ptr, std::array<std::uint32_t, dim> const& nCells)
        : ptr_{ptr}
        , nCells_{nCells}
    {
    }

    explicit NdArrayView(std::vector<std::decay_t<DataType>> const& v,
                         std::array<std::uint32_t, dim> const& nbCell)
        : NdArrayView{v.data(), nbCell}
    {
    }

    template<typename... Indexes>
    NO_DISCARD DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes...));
    }

    template<typename... Indexes>
    NO_DISCARD DataType const& operator()(Indexes... indexes) const
    {
        return viewer::at(ptr_, nCells_, indexes...);
    }

    template<typename Index>
    NO_DISCARD DataType const& operator()(std::array<Index, dim> const& indexes) const
    {
        return viewer::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    NO_DISCARD DataType& operator()(std::array<Index, dim> const& indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes));
    }

    NO_DISCARD auto data() const { return ptr_; }
    NO_DISCARD auto size() const
    {
        return std::accumulate(nCells_.begin(), nCells_.end(), std::size_t{1},
                               std::multiplies<std::size_t>());
    }
    NO_DISCARD auto shape() const { return nCells_; }

    void fill_from(NdArrayView const& that)
    {
        if (for_N_any<dim>([&](auto i) { return shape()[i] != that.shape()[i]; }))
            throw std::runtime_error("ArrayView::fill_from: Incompatible input shape");
        std::copy(that.data(), that.data() + size(), data());
    }

    NO_DISCARD auto begin() const { return ptr_; }
    NO_DISCARD auto begin() { return ptr_; }

    NO_DISCARD auto end() const { return ptr_ + size(); }
    NO_DISCARD auto end() { return ptr_ + size(); }

    void zero() { fill(0); }
    auto& fill(DataType const& v)
    {
        std::fill(begin(), end(), v);
        return *this;
    }

    void setBuffer(pointer_type ptr) { ptr_ = ptr; }
    void setShape(std::array<std::uint32_t, dim> const nCells) { nCells_ = nCells; }

private:
    pointer_type ptr_ = nullptr;
    std::array<std::uint32_t, dim> nCells_;
};


template<bool c_ordering = true, typename DataType, std::size_t dim>
auto make_array_view(DataType* data, std::array<std::uint32_t, dim> const& shape)
{
    return NdArrayView<dim, DataType, c_ordering>{data, shape};
}

template<bool c_ordering = true, typename DataType, std::size_t dim>
auto make_array_view(DataType const* const data, std::array<std::uint32_t, dim> shape)
{
    return NdArrayView<dim, DataType const, c_ordering>{data, shape};
}



template<std::size_t dim, typename DataType = double, bool c_ordering = true>
class NdArrayVector
{
public:
    static constexpr bool is_contiguous = 1;
    static std::size_t const dimension  = dim;
    using type                          = DataType;

    NdArrayVector() = delete;

    template<FloatingPoint U = DataType>
    explicit NdArrayVector(std::array<std::uint32_t, dim> const& ncells,
                           type const& value = static_cast<U>(std::nan("")))
        : nCells_{ncells}
        , data_(std::accumulate(ncells.begin(), ncells.end(), 1, std::multiplies<std::size_t>()),
                value)
    {
    }

    template<FloatingPoint U = DataType, typename... Nodes>
    explicit NdArrayVector(Nodes... nodes)
        requires(std::is_integral_v<Nodes> && ...)
        : NdArrayVector{std::array{nodes...}}
    {
        static_assert(sizeof...(Nodes) == dim);
    }

    explicit NdArrayVector(std::array<std::uint32_t, dim> const& ncells)
        requires(!FloatingPoint<DataType>)
        : nCells_{ncells}
        , data_(std::accumulate(ncells.begin(), ncells.end(), 1, std::multiplies<std::size_t>()))
    {
    }


    template<typename... Nodes>
        requires(!FloatingPoint<DataType>)
    explicit NdArrayVector(Nodes... nodes)
        requires(std::is_integral_v<Nodes> && ...)
        : NdArrayVector{std::array{nodes...}}
    {
        static_assert(sizeof...(Nodes) == dim);
    }


    NdArrayVector(NdArrayVector const& source)            = default;
    NdArrayVector(NdArrayVector&& source)                 = default;
    NdArrayVector& operator=(NdArrayVector const& source) = default;
    NdArrayVector& operator=(NdArrayVector&& source)      = default;

    NO_DISCARD auto data() const { return data_.data(); }
    NO_DISCARD auto data() { return data_.data(); }

    NO_DISCARD auto size() const { return data_.size(); }

    NO_DISCARD auto begin() const { return std::begin(data_); }
    NO_DISCARD auto begin() { return std::begin(data_); }

    NO_DISCARD auto end() const { return std::end(data_); }
    NO_DISCARD auto end() { return std::end(data_); }


    template<typename... Indexes>
    NO_DISCARD DataType const& operator()(Indexes... indexes) const
    {
        return NdArrayViewer<dim, c_ordering>::at(data_.data(), nCells_, indexes...);
    }

    template<typename... Indexes>
    NO_DISCARD DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes...));
    }

    template<typename Index>
    NO_DISCARD DataType const& operator()(std::array<Index, dim> const& indexes) const
    {
        return NdArrayViewer<dim, c_ordering>::at(data_.data(), nCells_, indexes);
    }

    template<typename Index>
    NO_DISCARD DataType& operator()(std::array<Index, dim> const& indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes));
    }


    NO_DISCARD auto& shape() const { return nCells_; }

    template<typename Mask>
    NO_DISCARD auto operator[](Mask&& mask)
    {
        return MaskedView{*this, std::forward<Mask>(mask)};
    }


    NO_DISCARD auto& vector() { return data_; }
    NO_DISCARD auto& vector() const { return data_; }


private:
    std::array<std::uint32_t, dim> nCells_;
    std::vector<DataType> data_;
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

        for (auto i = min_; i <= max_; ++i)
            array(i) = val;

        for (auto i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            array(i) = val;
    }

    template<typename Array>
    void fill2D(Array& array, typename Array::type val) const
    {
        auto shape = array.shape();

        // left border
        for (auto i = min_; i <= max_; ++i)
            for (auto j = min_; j <= shape[1] - 1 - max_; ++j)
                array(i, j) = val;

        // right border
        for (auto i = shape[0] - 1 - max_; i <= shape[0] - 1 - min_; ++i)
            for (auto j = min_; j <= shape[1] - 1 - max_; ++j)
                array(i, j) = val;


        for (auto i = min_; i <= shape[0] - 1 - min_; ++i)
        {
            // bottom border
            for (auto j = min_; j <= max_; ++j)
                array(i, j) = val;

            // top border
            for (auto j = shape[1] - 1 - max_; j <= shape[1] - 1 - min_; ++j)
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
    std::size_t min_, max_;
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
            outer(outer.xstart(), outer.ystart(), k) = inner(outer.xstart(), inner.ystart(), k);
            outer(outer.xstart(), outer.yend(), k)   = inner(outer.xstart(), inner.yend(), k);

            outer(outer.xend(), outer.ystart(), k) = inner(outer.xend(), inner.ystart(), k);
            outer(outer.xend(), outer.yend(), k)   = inner(outer.xend(), inner.yend(), k);
        }

        for (auto j = inner.ystart(); j <= inner.yend(); ++j)
        {
            outer(outer.xstart(), j, outer.zstart()) = inner(outer.xstart(), j, inner.zstart());
            outer(outer.xstart(), j, outer.zend())   = inner(outer.xstart(), j, inner.zend());

            outer(outer.xend(), j, outer.zstart()) = inner(outer.xend(), j, inner.zstart());
            outer(outer.xend(), j, outer.zend())   = inner(outer.xend(), j, inner.zend());
        }

        for (auto i = inner.xstart(); i <= inner.xend(); ++i)
        {
            outer(i, outer.ystart(), outer.zstart()) = inner(i, outer.ystart(), inner.zstart());
            outer(i, outer.ystart(), outer.zend())   = inner(i, outer.ystart(), inner.zend());

            outer(i, outer.yend(), outer.zstart()) = inner(i, outer.yend(), inner.zstart());
            outer(i, outer.yend(), outer.zend())   = inner(i, outer.yend(), inner.zend());
        }

        auto corner = [&](auto xouter, auto xinner) {
            outer(xouter, outer.ystart(), outer.zstart())
                = inner(xinner, outer.ystart(), inner.zstart());

            outer(xouter, outer.ystart(), outer.zend())
                = inner(xinner, outer.ystart(), inner.zend());

            outer(xouter, outer.yend(), outer.zstart())
                = inner(xinner, outer.yend(), inner.zstart());
            outer(xouter, outer.yend(), outer.zend()) = inner(xinner, outer.yend(), inner.zend());
        };

        corner(outer.xstart(), inner.xstart());
        corner(outer.xend(), inner.xend());
    }
}

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_HPP
