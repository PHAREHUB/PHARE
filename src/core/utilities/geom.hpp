#ifndef PHARE_CORE_UTILITIES_GEOM_HPP
#define PHARE_CORE_UTILITIES_GEOM_HPP

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

namespace PHARE::core
{
template<typename B, typename P, std::size_t dim = B::dimension>
std::int8_t normalized_icell_relative_to_closest_edge(B const& box, P const& point)
{
    // >= 0 = in domain     / 0 = first domain cell
    // <  0 = not in domain
    auto m = point;
    for (std::uint16_t i = 0; i < dim; ++i)
        m[i] = (box.lower[i] / 2) + (box.upper[i] / 2);
    for (std::uint16_t i = 0; i < dim; ++i)
    {
        if (point[i] < m[i])
            m[i] = point[i] - (box.lower[i]);
        else if (point[i] > m[i])
            m[i] = box.upper[i] - (point[i]);
    }

    std::int16_t min = *std::min_element(m.data(), m.data() + dim);
    if (min < std::numeric_limits<int8_t>::min() or min > std::numeric_limits<int8_t>::max())
        throw std::runtime_error("no");
    return min;
}


template<typename Box_t, std::size_t impl_ = 0>
struct PointInBoxChecker
{
    template<typename Point_t>
    bool operator()(Point_t const& p);

    Box_t const& box;
};

template<typename Box_t, typename Point_t>
bool point_in_box_0(Box_t const& box, Point_t const& point)
{
    auto isIn1D = [](auto const& pos, auto const& lower, auto const& upper) {
        return pos >= lower && pos <= upper;
    };
    bool pointInBox = true;
    for (auto iDim = 0u; iDim < point.size(); ++iDim)
        pointInBox = pointInBox && isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
    return pointInBox;
}
template<typename Box_t, typename Point_t>
bool point_in_box_1(Box_t const& box, Point_t const& point)
{
    auto isIn1D = [](auto const& pos, auto const& lower, auto const& upper) {
        return pos >= lower & pos <= upper;
    };
    bool pointInBox = true;
    for (auto iDim = 0u; iDim < point.size(); ++iDim)
        pointInBox &= isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
    return pointInBox;
}
template<typename Box_t, typename Point_t>
bool point_in_box_2(Box_t const& box, Point_t const& point)
{
    for (auto iDim = 0u; iDim < point.size(); ++iDim)
        if (!(point[iDim] < box.lower[iDim] | point[iDim] > box.upper[iDim]))
            return true;
    return false;
}


auto static inline greater(int const a, int const& b)
{
    return (a - ((a - b) & ((a - b) >> (sizeof(int) * CHAR_BIT - 1))));
}
auto static inline is_greater(int const& a, int const& b)
{
    return ((a != b) & (greater(a, b) == a));
}
auto static inline is_or_greater(int const& a, int const& b)
{
    return greater(a, b) == a;
}

template<typename Box_t, typename Point_t>
bool point_in_box_3(Box_t const& box, Point_t const& point)
{
    auto isIn1D = [](auto const& pos, auto const& lower, auto const& upper) {
        return is_or_greater(pos, lower) & is_greater(upper, pos);
    };
    bool pointInBox = true;
    for (auto iDim = 0u; iDim < point.size(); ++iDim)
        pointInBox &= isIn1D(point[iDim], box.lower[iDim], box.upper[iDim]);
    return pointInBox;
}

template<typename Box_t, std::size_t impl>
template<typename Point_t>
bool PointInBoxChecker<Box_t, impl>::operator()(Point_t const& p)
{
    static_assert(impl >= 0 and impl < 4);

    if constexpr (impl == 0)
        return point_in_box_0(box, p);
    if constexpr (impl == 1)
        return point_in_box_1(box, p);
    if constexpr (impl == 2)
        return point_in_box_2(box, p);
    if constexpr (impl == 3)
        return point_in_box_3(box, p);
}


} // namespace PHARE::core

#endif /* PHARE_CORE_UTILITIES_GEOM_HPP */
