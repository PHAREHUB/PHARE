
#ifndef PHARE_CORE_UTILITIES_BOX_BOX_H
#define PHARE_CORE_UTILITIES_BOX_BOX_H

#include <cstddef>

#include "core/utilities/point/point.h"
#include "core/utilities/meta/meta_utilities.h"
#include "core/def.h"

namespace PHARE
{
namespace core
{
    /** Represents a 1D, 2D or 3D box of integer or floating point
     * points.
     */
    template<typename Type, std::size_t dim>
    struct Box
    {
        Point<Type, dim> lower;
        Point<Type, dim> upper;

        Box() = default;

        template<typename T, std::size_t s>
        Box(Point<T, s> _lower, Point<T, s> _upper)
            : lower{_lower}
            , upper{_upper}
        {
        }

        bool operator==(Box const& box) const { return box.lower == lower && box.upper == upper; }

        bool isEmpty() const { return (*this) == Box{}; }

        auto nbrItems(std::size_t dir) const { return upper[dir] - lower[dir]; }


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

        using type = Type;
    };

    template<typename T, std::size_t s>
    Box(Point<T, s> lower, Point<T, s> upper) -> Box<T, s>;


    template<typename T, std::size_t dim>
    bool sameSize(Box<T, dim> const& box1, Box<T, dim> const& box2)
    {
        static_assert(std::is_integral_v<T>,
                      "this function is only valid for integral type of Point");
        bool same = true;
        for (auto i = 0u; i < dim; ++i)
        {
            same &= ((box1.upper[i] - box1.lower[i]) == (box2.upper[i] - box2.lower[i]));
        }
        return same;
    }



    /** this overload of isIn takes a Point and a Container of boxes
     * and returns true if the Point is at least in one of the boxes.
     * Returns occurs at the first box the point is in.
     */
    template<typename Point, typename BoxContainer, is_iterable<BoxContainer> = dummy::value>
    bool isIn(Point const& point, BoxContainer const& boxes) _PHARE_FN_SIG_
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
    bool isIn(Point const& point,
              Box<typename Point::type, Point::dimension> const& box) _PHARE_FN_SIG_
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

} // namespace core
} // namespace PHARE

#endif
