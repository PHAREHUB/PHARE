#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SPANS_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SPANS_HPP

#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace PHARE::core
{

/*  Span by span walk of a box.

    A span is the run of cells along the last (fastest varying) direction, the
    other indices fixed, given as (start point, size). It is an index range, not
    a Span<T> (core/utilities/span.hpp): it holds no memory. In a C ordered array
    (NdArrayVector) the cells of a span are contiguous in memory: a caller computes
    one array offset per span, then runs a plain loop over it, which the compiler
    vectorizes.

        for_each_span(box, [&](auto const& start, auto const size) { ... });
        for_each_span(dst_box, src_box, [&](auto const& dst_start, auto const& src_start,
                                            auto const size) { ... });
        for_each_point(box, [&](auto const& point) { ... });

    The walk is plain nested loops: going to the next span is just ++y,
    with no iterator object and no y -> x carry test.
*/


/** Walks `box` span by span, in memory order: calls fn(start, size) where start is
 * the first point of the span and size its number of cells.
 * `start` is only valid for the duration of the call.
 * A box with lower > upper in any direction (e.g. an empty intersection) has no spans.
 */
template<typename T, std::size_t dim>
void for_each_span(Box<T, dim> const& box, auto&& fn)
{
    static_assert(dim >= 1 and dim <= 3);

    // one past the last cell in each direction. `x < end` rather than `x <= upper`
    // also makes an unsigned box with upper = lower - 1 = UINT_MAX empty
    auto const end = box.upper + 1;
    auto const& lo = box.lower;

    if (end[dim - 1] <= lo[dim - 1])
        return; // spans are empty. An empty outer direction makes its own loop run zero times

    std::uint32_t const size = end[dim - 1] - lo[dim - 1];

    if constexpr (dim == 1)
        fn(lo, size);

    if constexpr (dim == 2)
        for (T x = lo[0]; x < end[0]; ++x)
            fn(Point<T, 2>{x, lo[1]}, size);

    if constexpr (dim == 3)
        for (T x = lo[0]; x < end[0]; ++x)
            for (T y = lo[1]; y < end[1]; ++y)
                fn(Point<T, 3>{x, y, lo[2]}, size);
}


/** Walks two same-shaped boxes span by span, in lockstep and in memory order:
 * calls fn(a_start, b_start, size) where a_start and b_start are the first
 * points of matching spans of `a` and `b`, and size is their number of cells.
 * The start points are only valid for the duration of the call.
 */
template<typename T, std::size_t dim>
void for_each_span(Box<T, dim> const& a, Box<T, dim> const& b, auto&& fn)
{
    assert(a.shape() == b.shape()); // b is only read through its lower corner

    auto const shift = b.lower - a.lower;
    for_each_span(a, [&](auto const& start, auto const size) { fn(start, start + shift, size); });
}


/** Calls fn(point) for every point of `box`, in memory order.
 * `point` is only valid for the duration of the call.
 */
template<typename T, std::size_t dim>
void for_each_point(Box<T, dim> const& box, auto&& fn)
{
    for_each_span(box, [&](auto const& start, std::uint32_t const size) {
        auto point = start;
        // `<` and not `!=`: the compiler then knows the 32 bit index cannot wrap,
        // so array offsets stay affine in i and the loop can be vectorized
        for (T i = start[dim - 1], end = i + size; i < end; ++i)
        {
            point[dim - 1] = i;
            fn(point);
        }
    });
}


} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_BOX_BOX_SPANS_HPP
