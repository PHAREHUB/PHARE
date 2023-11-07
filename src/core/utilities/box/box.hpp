#ifndef PHARE_CORE_UTILITIES_BOX_BOX_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/def.hpp"

#include <cstddef>
#include <algorithm>
#include <optional>
#include <iostream>

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
    Point<Type, dim> lower;
    Point<Type, dim> upper;

    Box() = default;

    constexpr Box(std::array<Type, dim> _lower, std::array<Type, dim> _upper)
        : lower{_lower}
        , upper{_upper}
    {
    }

    template<typename T, std::size_t s>
    Box(Point<T, s> _lower, Point<T, s> _upper)
        : lower{_lower}
        , upper{_upper}
    {
    }

    NO_DISCARD bool operator==(Box const& box) const
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
        for (auto& c : lower)
        {
            c -= size;
        }
        for (auto& c : upper)
        {
            c += size;
        }
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
    Point<Type, dim> operator*() { return index_; }


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
    box_iterator& operator++() const
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
template<typename Point>
NO_DISCARD bool isIn(Point const& point,
                     Box<typename Point::value_type, Point::dimension> const& box)
{
    auto isIn1D = [](typename Point::value_type pos, typename Point::value_type lower,
                     typename Point::value_type upper) { return pos >= lower && pos <= upper; };

    bool pointInBox = true;

    for (auto iDim = 0u; iDim < Point::dimension; ++iDim)
    {
        pointInBox = pointInBox && isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
    }
    if (pointInBox)
        return pointInBox;

    return false;
}


template<typename Type, std::size_t dim>
Box<Type, dim> grow(Box<Type, dim> const& box, Type const& size)
{
    auto copy{box};
    copy.grow(size);
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



} // namespace PHARE::core

#endif
