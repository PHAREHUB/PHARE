#ifndef PHARE_CORE_GRID_GRIDLAYOUTDEFS_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTDEFS_HPP

#include "core/utilities/point/point.hpp"

#include <cstddef>

namespace PHARE
{
namespace core
{
    enum class Direction { X, Y, Z };


    enum class QtyCentering : std::uint16_t { primal = 0, dual = 1 };


    template<std::size_t dim>
    struct WeightPoint
    {
        Point<int, dim> indexes{};
        double coef{0.0};
    };

    // using LinearCombination = std::vector<WeightPoint>;

    enum class Layout { Yee };


} // namespace core
} // namespace PHARE

#endif // GRIDLAYOUTDEFS_HPP
