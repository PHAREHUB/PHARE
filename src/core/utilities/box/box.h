
#ifndef PHARE_CORE_UTILITIES_BOX_BOX_H
#define PHARE_CORE_UTILITIES_BOX_BOX_H

#include <cstddef>
#include <algorithm>

#include "core/utilities/point/point.h"
#include "core/utilities/meta/meta_utilities.h"


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
    std::array<Type, dim> shape;

    Box() = default;

    Box(std::array<Type, dim> _lower, std::array<Type, dim> _upper)
        : lower{_lower}
        , upper{_upper}
    {
        make_shape_();
    }

    template<typename T, std::size_t s>
    Box(Point<T, s> _lower, Point<T, s> _upper)
        : lower{_lower}
        , upper{_upper}
    {
        make_shape_();
    }

    bool operator==(Box const& box) const { return box.lower == lower && box.upper == upper; }

    auto operator*(Box const& other) const
    {
        Box intersection{other};
        for (auto idim = 0u; idim < dim; ++idim)
        {
            intersection.lower[idim] = std::max(lower[idim], other.lower[idim]);
            intersection.upper[idim] = std::min(upper[idim], other.upper[idim]);
            assert(intersection.lower[idim] <= intersection.upper[idim]);
        }
        return intersection;
    }

    bool isEmpty() const { return (*this) == Box{}; }

    auto nbrItems(std::size_t dir) const { return shape[dir]; }

    auto size() const
    {
        return std::accumulate(std::begin(shape), std::end(shape), 1, std::multiplies<Type>{});
    }

    void grow(Type size)
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

    using iterator = box_iterator<Type, dim>;
    auto begin() { return iterator{this, lower}; }

    //   // since the 1D scan of the multidimensional box is done assuming C ordering
    //   // the end (in the sense of container.end()) is one beyond last for the last
    //   // direction only, previous dimensions have not reached the end.
    auto begin() const { return iterator{this, lower}; }
    auto end()
    {
        static_assert(dim <= 3 and dim > 0);
        // following could maybe be a one liner?
        if constexpr (dim == 1)
        {
            return iterator{this, {upper[0] + 1}};
        }
        else if constexpr (dim == 2)
        {
            return iterator{this, {upper[0], upper[1] + 1}};
        }
        else
        {
            return iterator{this, {upper[0], upper[1], upper[2] + 1}};
        }
    }

    auto end() const
    {
        static_assert(dim <= 3 and dim > 0);
        if constexpr (dim == 1)
        {
            return iterator{this, {upper[0] + 1}};
        }
        else if constexpr (dim == 2)
        {
            return iterator{this, {upper[0], upper[1] + 1}};
        }
        else
        {
            return iterator{this, {upper[0], upper[1], upper[2] + 1}};
        }
    }
    using type = Type;

private:
    void make_shape_()
    { //
        std::transform(std::begin(lower), std::end(lower), std::begin(upper), std::begin(shape),
                       [](auto const& lo, auto const& up) { return up - lo + 1; });
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

    box_iterator operator++()
    {
        // this ++ operator assumes C ordering
        // this means that incrementing
        // cell at (i,j,k) will go to cell
        // (i,j,k+1).
        // so the first thing we need to do
        // is incrementing the last dimension index
        index_[dim - 1]++;

        // then it is possible that this increment has
        // made the index k reaching last index+1
        // if it is the case, we need to reset it to
        // the lower cell of the box in that direction
        // and increment j. This is only done if j itself
        // has not reached the last index
        // if j reaches last+1, then set it to lower j
        // and increment i.
        for (auto idim = dim - 1; idim > 0; idim--)
        {
            // note : be sure this works in 3D too...
            if (index_[idim] == box_->upper[idim] + 1 and index_[idim - 1] < box_->upper[idim - 1])
            {
                index_[idim] = box_->lower[idim];
                index_[idim - 1]++;
            }
        }
        return *this;
    }
    box_iterator operator++() const
    {
        index_[dim - 1]++;
        for (auto idim = dim - 1; idim > 0; idim--)
        {
            if (index_[idim] == box_->upper[idim] + 1 and index_[idim - 1] < box_->upper[idim - 1])
            {
                index_[idim] = box_->lower[idim];
                index_[idim - 1]++;
            }
        }
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


    static_assert(
        std::is_same<typename Point::type, typename BoxContainer::value_type::type>::value,
        "Box and Point should have the same data type");


    auto isIn1D = [](typename Point::type pos, typename Point::type lower,
                     typename Point::type upper) { return pos >= lower && pos <= upper; };

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
bool isIn(Point const& point, Box<typename Point::type, Point::dimension> const& box)
{
    auto isIn1D = [](typename Point::type pos, typename Point::type lower,
                     typename Point::type upper) { return pos >= lower && pos <= upper; };

    bool pointInBox = true;

    for (auto iDim = 0u; iDim < Point::dimension; ++iDim)
    {
        pointInBox = pointInBox && isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
    }
    if (pointInBox)
        return pointInBox;

    return false;
}

} // namespace PHARE::core

#endif
