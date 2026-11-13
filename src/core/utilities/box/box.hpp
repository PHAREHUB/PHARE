#ifndef PHARE_CORE_UTILITIES_BOX_BOX_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_HPP

#include "core/def.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"

#include <cstddef>
#include <optional>
#include <iostream>
#include <algorithm>

namespace PHARE::core
{
/** Represents a 1D, 2D or 3D box of integer or floating point
 * points.
 */
template<typename Type, std::size_t dim>
struct Box
{
    using value_type                = Type;
    static constexpr auto dimension = dim;


    Point<Type, dim> lower;
    Point<Type, dim> upper;

    Box() = default;

    constexpr Box(std::array<Type, dim> const& _lower, std::array<Type, dim> const& _upper)
        : lower{_lower}
        , upper{_upper}
    {
    }

    template<typename T, std::size_t s>
    Box(Point<T, s> const& _lower, Point<T, s> const& _upper)
        : lower{_lower}
        , upper{_upper}
    {
    }

    template<typename T2>
    NO_DISCARD bool operator==(Box<T2, dim> const& box) const
    {
        return box.lower == lower && box.upper == upper;
    }

    NO_DISCARD auto operator*(Box const& other) const
    {
        Box intersection{other};
        for (auto idim = 0u; idim < dim; ++idim)
        {
            intersection.lower[idim] = std::max(lower[idim], other.lower[idim]);
            intersection.upper[idim] = std::min(upper[idim], other.upper[idim]);
            if (intersection.lower[idim] > intersection.upper[idim])
                return std::optional<Box>{std::nullopt};
        }
        return std::optional<Box>{intersection};
    }

    NO_DISCARD auto operator-(Box const& that) const
    {
        return Box{lower - that.lower, upper - that.upper};
    }


    template<typename Modifier>
    auto& grow(Modifier const& size)
    {
        lower -= size;
        upper += size;
        return *this;
    }

    template<typename Modifier>
    auto& shrink(Modifier const& size)
    {
        lower += size;
        upper -= size;
        return *this;
    }


    NO_DISCARD auto shape() const { return upper - lower + 1; }
    NO_DISCARD auto size() const { return core::product(shape(), std::size_t{1}); }


    NO_DISCARD constexpr static std::size_t nbrRemainBoxes()
    {
        if constexpr (dim == 1)
        {
            return 2;
        }
        else if constexpr (dim == 2)
        {
            return 4;
        }
        else
            return 6;
    }
};


template<typename T, std::size_t s>
Box(Point<T, s> lower, Point<T, s> upper) -> Box<T, s>;




template<template<typename, std::size_t> typename Point, typename Type, std::size_t SIZE>
NO_DISCARD bool isIn(Point<Type, SIZE> const& point, Box<Type, SIZE> const& box)
{
    for (std::size_t iDim = 0; iDim < SIZE; ++iDim)
        if (point[iDim] < box.lower[iDim] || point[iDim] > box.upper[iDim])
            return false;
    return true;
}

template<Spannable Boxes>
NO_DISCARD bool isIn(auto const& point, Boxes const& boxes)
{
    for (auto const& box : boxes)
        if (isIn(point, box))
            return true;
    return false;
}



template<typename Type, std::size_t dim, typename OType>
Box<Type, dim> grow(Box<Type, dim> const& box, OType const& size)
{
    auto copy{box};
    copy.grow(size);
    return copy;
}

template<typename Type, std::size_t dim, typename T2>
NO_DISCARD Box<Type, dim> shrink(Box<Type, dim> const& box, T2 const& size)
{
    auto copy{box};
    copy.shrink(size);
    return copy;
}


template<typename Type, std::size_t dim, typename Shifter>
NO_DISCARD Box<Type, dim> shift(Box<Type, dim> const& box, Shifter const& offset)
{
    auto copy{box};
    copy.lower += offset;
    copy.upper += offset;
    return copy;
}



template<typename Type, std::size_t dim>
auto& operator<<(std::ostream& os, Box<Type, dim> const& box)
{
    os << "Box<Type," << dim << "> : ( ";
    for (auto& c : box.lower)
        os << c << " ";
    os << ")-->( ";
    for (auto& c : box.upper)
        os << c << " ";
    os << ")";
    return os;
}


} // namespace PHARE::core

#endif
