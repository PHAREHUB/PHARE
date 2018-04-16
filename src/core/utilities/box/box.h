
#ifndef PHARE_CORE_UTILITIES_BOX_BOX_H
#define PHARE_CORE_UTILITIES_BOX_BOX_H

#include <cstddef>

#include "utilities/point/point.h"


namespace PHARE
{
template<typename Type, std::size_t dim>
struct Box
{
    Point<Type, dim> lower;
    Point<Type, dim> upper;

    using type = Type;
};


template<typename Point, typename BoxContainer>
bool isIn(Point const& point, BoxContainer const& boxes)
{
    if (boxes.size() == 0)
        return false;


    static_assert(
        std::is_same<typename Point::type, typename BoxContainer::value_type::type>::value,
        "Box and Point should have the same data type");


    auto isIn1D = [](typename Point::type pos, typename Point::type lower,
                     typename Point::type upper) { return pos >= lower && pos < upper; };

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


} // namespace PHARE

#endif
