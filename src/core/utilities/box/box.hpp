#ifndef PHARE_CORE_UTILITIES_BOX_BOX_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include <limits>
#include <cstddef>
#include <optional>
#include <iostream>
#include <algorithm>

namespace PHARE::core
{
template<typename Type, std::size_t dim>
class box_iterator;



/** Represents a 1D, 2D or 3D box of integer or floating point
 * points.
 */
template<typename Type, std::size_t dim>
struct Box
{
    static size_t const dimension = dim;
    using Point_t                 = Point<Type, dim>;


    Point_t lower;
    Point_t upper;

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
        assert(lower <= upper);
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


    NO_DISCARD bool isEmpty() const { return (*this) == Box{}; }

    void grow(Type const& size)
    {
        assert(size >= 0);
        lower -= size;
        upper += size;
    }


    template<typename Size>
    auto& grow(std::array<Size, dim> const& size)
    {
        lower -= size;
        upper += size;
        return *this;
    }

    auto& shrink(Type const& size)
    {
        assert(size >= 0);
        lower += size;
        upper -= size;
        return *this;
    }

    template<typename Size>
    auto& shrink(std::array<Size, dim> const& size)
    {
        lower += size;
        upper -= size;
        return *this;
    }


    NO_DISCARD auto shape() const { return upper - lower + 1; }
    NO_DISCARD auto size() const { return core::product(shape()); }


    using iterator = box_iterator<Type, dim>;
    NO_DISCARD auto begin() { return iterator{this, lower}; }

    //   // since the 1D scan of the multidimensional box is done assuming C ordering
    //   // the end (in the sense of container.end()) is one beyond last for the last
    //   // direction only, previous dimensions have not reached the end.
    NO_DISCARD auto begin() const { return iterator{this, lower}; }
    NO_DISCARD auto end()
    {
        static_assert(dim <= 3 and dim > 0);
        // following could maybe be a one liner?
        if constexpr (dim == 1)
        {
            return iterator{this, {upper[0] + 1}};
        }
        else if constexpr (dim == 2)
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1}};
        }
        else
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1, upper[2] + 1}};
        }
    }

    NO_DISCARD auto end() const
    {
        static_assert(dim <= 3 and dim > 0);
        if constexpr (dim == 1)
        {
            return iterator{this, {upper[0] + 1}};
        }
        else if constexpr (dim == 2)
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1}};
        }
        else
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1, upper[2] + 1}};
        }
    }
    using value_type = Type;


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


    std::vector<Box> remove(Box const& that) const;
    Box merge(Box const& that) const;
    Box remove_edge(Box const& that) const;
};

template<typename Type, std::size_t dim>
class box_iterator
{
public:
    box_iterator(Box<Type, dim> const* box, Point<Type, dim> index = Point<Type, dim>{})
        : box_{box}
        , index_{index}
    {
    }

public:
    auto& operator*() const { return index_; }
    auto operator->() const { return &index_; }

    void increment(std::size_t idim)
    {
        index_[idim]++;
        if (idim == 0)
            return;
        if (index_[idim] == box_->upper[idim] + 1)
        {
            increment(idim - 1);
            if (index_[idim - 1] <= box_->upper[idim - 1])
                index_[idim] = box_->lower[idim];
        }
    }

    box_iterator& operator++()
    {
        increment(dim - 1);
        return *this;
    }

    bool operator!=(box_iterator const& other) const
    {
        return box_ != other.box_ or index_ != other.index_;
    }


private:
    Box<Type, dim> const* box_;
    Point<Type, dim> index_;
};


template<typename T, std::size_t s>
Box(Point<T, s> lower, Point<T, s> upper) -> Box<T, s>;



/** this overload of isIn takes a Point and a Container of boxes
 * and returns true if the Point is at least in one of the boxes.
 * Returns occurs at the first box the point is in.
 */
template<typename Point, typename BoxContainer, is_iterable<BoxContainer> = dummy::value>
bool isIn(Point const& point, BoxContainer const& boxes)
{
    if (boxes.size() == 0)
        return false;


    static_assert(std::is_same<typename Point::value_type,
                               typename BoxContainer::value_type::value_type>::value,
                  "Box and Point should have the same data type");


    auto isIn1D = [](typename Point::value_type pos, typename Point::value_type lower,
                     typename Point::value_type upper) { return pos >= lower && pos <= upper; };

    for (auto const& box : boxes)
    {
        bool pointInBox = true;

        for (auto iDim = 0u; iDim < Point::dimension; ++iDim)
        {
            pointInBox = pointInBox && isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
        }
        if (pointInBox)
            return pointInBox;
    }

    return false;
}

/** This overload of isIn does the same as the one above but takes only
 * one box.
 */
template<template<typename, std::size_t> typename Point, typename Type, std::size_t SIZE>
NO_DISCARD bool isIn(Point<Type, SIZE> const& point, Box<Type, SIZE> const& box)
{
    auto isIn1D = [](auto const pos, auto const lower, auto const upper) {
        return pos >= lower && pos <= upper;
    };

    bool pointInBox = true;

    for (auto iDim = 0u; iDim < SIZE; ++iDim)
        pointInBox = pointInBox && isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
    if (pointInBox)
        return pointInBox;

    return false;
}

template<typename Particle, typename Type>
NO_DISCARD auto isIn(Particle const& particle, Box<Type, Particle::dimension> const& box)
    -> decltype(isIn(particle.iCell, box), bool())
{
    return isIn(particle.iCell, box);
}

template<typename... Args>
NO_DISCARD auto isIn(Args const&&... args)
{
    return isIn(args...);
}


template<template<typename, std::size_t> typename Point, typename Type, std::size_t SIZE>
struct InBox
{
    bool operator()(auto const& arg) const { return isIn(arg, box); }

    Box<Type, SIZE> const& box;
};


template<typename Type, std::size_t dim, typename OType>
NO_DISCARD Box<Type, dim> grow(Box<Type, dim> const& box, OType const& size)
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

template<typename Type, std::size_t dim>
NO_DISCARD Box<Type, dim> emptyBox()
{
    return Box<Type, dim>{};
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



template<typename Type, std::size_t dim>
std::vector<Box<Type, dim>> Box<Type, dim>::remove(Box<Type, dim> const& to_remove) const
{
    using box_t = Box<Type, dim>;
    using _m    = std::unordered_map<std::uint16_t, std::uint32_t>;

    auto const box = *this;

    auto overlap = box * to_remove;

    if (not overlap)
        return std::vector{box};

    auto copy = [](auto cpy, auto const& replace) {
        for (auto const& [i, v] : replace)
            cpy[i] = v;
        return cpy;
    };

    auto intersection = *overlap;

    std::unordered_map<std::string, box_t> boxes;

    if (intersection.lower[0] > box.lower[0])
        boxes["left"] = Box(box.lower, copy(box.upper, _m{{0, intersection.lower[0] - 1}}));
    if (intersection.upper[0] < box.upper[0])
        boxes["right"] = box_t{copy(box.lower, _m{{0, intersection.upper[0] + 1}}), box.upper};

    [[maybe_unused]] Type minx = 0, maxx = 0;
    if constexpr (dim > 1)
    {
        minx = boxes.count("left") > 0 ? intersection.lower[0] : box.lower[0];
        maxx = boxes.count("right") > 0 ? intersection.upper[0] : box.upper[0];

        if (intersection.lower[1] > box.lower[1])
            boxes["down"] = box_t{copy(box.lower, _m{{0, minx}}),
                                  copy(box.upper, _m{{0, maxx}, {1, intersection.lower[1] - 1}})};

        if (intersection.upper[1] < box.upper[1])
            boxes["up"] = Box(copy(box.lower, _m{{0, minx}, {1, intersection.upper[1] + 1}}),
                              copy(box.upper, _m{{0, maxx}}));
    }

    if constexpr (dim > 2)
    {
        Type miny = boxes.count("down") > 0 ? intersection.lower[1] : box.lower[1];
        Type maxy = boxes.count("up") > 0 ? intersection.upper[1] : box.upper[1];

        if (intersection.lower[2] > box.lower[2])
            boxes["back"] = Box(copy(box.lower, _m{{0, minx}, {1, miny}}),
                                copy(intersection.lower - 1, _m{{0, maxx}, {1, maxy}}));
        if (intersection.upper[2] < box.upper[2])
            boxes["front"] = Box(copy(intersection.upper + 1, _m{{0, minx}, {1, miny}}),
                                 copy(box.upper, _m{{0, maxx}, {1, maxy}}));
    }

    std::vector<box_t> remaining;
    for (auto const& [key, val] : boxes)
        remaining.emplace_back(val);
    return remaining;
}


template<typename Type, std::size_t dim>
Box<Type, dim> Box<Type, dim>::merge(Box<Type, dim> const& to_merge) const
{
    Box<Type, dim> merged{*this};
    for (std::size_t i = 0; i < dim; ++i)
    {
        merged.lower[i] = std::min(merged.lower[i], to_merge.lower[i]);
        merged.upper[i] = std::max(merged.upper[i], to_merge.upper[i]);
    }
    return merged;
}

template<typename Type, std::size_t dim>
Box<Type, dim> Box<Type, dim>::remove_edge(Box<Type, dim> const& to_remove) const
{
    auto const remaining = this->remove(to_remove);

    if (remaining.size() == 0)
        return *this;
    if (remaining.size() == 1)
        return remaining[0];

    Box<Type, dim> leftover{ConstArray<Type, dim>(std::numeric_limits<Type>::max()),
                            ConstArray<Type, dim>(std::numeric_limits<Type>::min())};

    for (auto const& rem : remaining)
        leftover = leftover.merge(rem);

    return leftover;
}

template<typename BoxHavers, typename Accessor>
bool any_overlaps_in(BoxHavers const& havers, Accessor&& fn)
{
    for (std::size_t i = 0; i < havers.size() - 1; ++i)
        for (std::size_t j = i + 1; j < havers.size(); ++j)
            if (auto overlap = fn(havers[i]) * fn(havers[j]))
                return true;
    return false;
}

template<typename Type, std::size_t dim>
NO_DISCARD Box<Type, dim> shift(Box<Type, dim> const& box, Type const& offset)
{
    auto copy{box};
    copy.lower += offset;
    copy.upper += offset;
    return copy;
}

template<template<typename, std::size_t> typename Point_t, typename Type, std::size_t dim>
NO_DISCARD Box<Type, dim> shift(Box<Type, dim> const& box, Point_t<Type, dim> const& offset)
{
    auto copy{box};
    for (std::uint8_t i = 0; i < dim; ++i)
        copy.lower[i] += offset[i], copy.upper[i] += offset[i];
    return copy;
}

template<std::uint8_t idx, typename Type, std::size_t dim>
NO_DISCARD Box<Type, dim> shift_idx(Box<Type, dim> const& box, Type const& offset)
{
    auto copy{box};
    copy.lower[idx] += offset;
    copy.upper[idx] += offset;
    return copy;
}

} // namespace PHARE::core

#endif
