#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SLABS_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SLABS_HPP

#include "core/utilities/box/box.hpp"

#include <cstddef>
#include <algorithm>
#include <cstdint>


namespace PHARE::core
{

/*
Span-wise iteration over an N dimensional Box.

Instead of visiting the box point by point, the box is decomposed into:
  - slabs: one per index of the slowest-varying dimension (3D only,
           1D/2D boxes are a single slab)
  - spans: within a slab, one per index of the second-fastest dimension;
           each span is the contiguous run of cells along the last
           (fastest-varying) dimension

    for (auto const& slab : make_box_slabs(box))   // 3D: one slab per x index
        for (auto const& [start, size] : slab)    // one span per y index
            ...  // span covers cells [start, start + size) along z

This gives callers whole contiguous spans to work with (memcpy, vectorizable
inner loops) rather than a callback per point.
*/


/** @brief One span of a box: the contiguous run of `size` cells from `start` along the
 * fastest-varying dimension.
 *
 * Binds as [start, size], or iterates its points: a point_iterator advances its own
 * copy of start, which unlike advancing start in place lets the compiler vectorize
 * loops over the span (2D points otherwise stay packed in one 64 bit register).
 */
template<typename T, std::size_t dim>
struct BoxSpan
{
    class point_iterator;

    point_iterator begin() const { return {start}; }
    T end() const { return start[dim - 1] + size; }

    Point<T, dim> start;
    std::uint32_t size;
};


/** @brief Entry point of the span-wise iteration: iterates the slabs of a box.
 * See make_box_slabs() for the usual way to construct one.
 */
template<typename T, std::size_t dim_>
class BoxSlabs
{
    auto constexpr static dim = dim_;

public:
    using value_type = T;

    class slab_iterator;
    class span_iterator;

    BoxSlabs(Box<T, dim> const& _box)
        : box{_box}
    {
    }

    slab_iterator begin() { return {box, slab_begin()}; }
    slab_iterator begin() const { return {box, slab_begin()}; }
    slab_iterator end() { return {box, slab_end()}; }
    slab_iterator end() const { return {box, slab_end()}; }

    T slab_begin() const
    {
        if constexpr (dim > 2)
            return box.lower[0];
        return 0;
    }
    T slab_end() const // clamped, see span_iterator::_span_size
    {
        if constexpr (dim > 2)
            return std::max<T>(box.lower[0], box.upper[0] + 1);
        return 1;
    }

    // number of slabs (not the number of cells in the box), 0 for an ill-formed box
    std::uint32_t size() const { return slab_end() - slab_begin(); }

protected:
    Box<T, dim> const box;
};


/** @brief Iterator over the spans of one slab.
 *
 * Dereferencing yields the current BoxSpan by value.
 *
 * Comparison operators only inspect span_idx, so comparing iterators made
 * from different boxes or slabs is meaningless.
 */
template<typename T, std::size_t dim_>
class BoxSlabs<T, dim_>::span_iterator
{
public:
    using value_type = T;

    span_iterator(Box<T, dim> const& box, T const slab_idx, T const span_idx)
        : box{box}
        , slab_idx{slab_idx}
        , span_idx{span_idx}
        , span_size{_span_size()}
    {
    }

    span_iterator& operator++()
    {
        ++span_idx;
        return *this;
    }
    span_iterator operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    span_iterator& operator--()
    {
        --span_idx;
        return *this;
    }
    span_iterator operator--(int) // postfix decrement
    {
        auto copy = *this;
        --(*this);
        return copy;
    }

    // starting point of the current span: fastest-varying coordinate is always
    // box.lower[dim - 1], since a span walks that whole dimension contiguously
    auto& point()
    {
        if constexpr (dim == 1)
            return (_point = {box.lower[0]});
        if constexpr (dim == 2)
            return (_point = {span_idx, box.lower[1]});
        if constexpr (dim == 3)
            return (_point = {slab_idx, span_idx, box.lower[2]});
    }

    BoxSpan<T, dim> operator*() { return {point(), span_size}; }

    bool operator==(span_iterator const& that) const { return span_idx == that.span_idx; }
    bool operator!=(span_iterator const& that) const { return span_idx != that.span_idx; }
    bool operator<(span_iterator const& that) const { return span_idx < that.span_idx; }
    bool operator>(span_iterator const& that) const { return span_idx > that.span_idx; }


protected:
    // ill-formed boxes (lower > upper, e.g. an empty SAMRAI intersection) are clamped to empty
    std::uint32_t _span_size() const
    {
        return box.upper[dim - 1] >= box.lower[dim - 1]
                   ? box.upper[dim - 1] - box.lower[dim - 1] + 1
                   : 0;
    }

    Box<T, dim> const& box;
    T const slab_idx;
    T span_idx;
    std::uint32_t const span_size;
    Point<T, dim> _point{};
};


/** @brief Iterator over the slabs of a box.
 *
 * A slab is a fixed slowest-varying index: dereferencing yields the slab_iterator
 * itself, whose begin()/end() are span_iterators over that slab's spans.
 */
template<typename T, std::size_t dim_>
class BoxSlabs<T, dim_>::slab_iterator
{
public:
    using value_type = T;

    slab_iterator(Box<T, dim> const& _box, T const _slab_idx)
        : box{_box}
        , slab_idx{_slab_idx}
    {
    }

    span_iterator begin() { return {box, slab_idx, span_begin()}; }
    span_iterator begin() const { return {box, slab_idx, span_begin()}; }
    span_iterator end() { return {box, slab_idx, span_end()}; }
    span_iterator end() const { return {box, slab_idx, span_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

    bool operator==(slab_iterator const& that) const { return slab_idx == that.slab_idx; }
    bool operator!=(slab_iterator const& that) const { return slab_idx != that.slab_idx; }
    bool operator<(slab_iterator const& that) const { return slab_idx < that.slab_idx; }
    bool operator>(slab_iterator const& that) const { return slab_idx > that.slab_idx; }

    slab_iterator& operator++()
    {
        ++slab_idx;
        return *this;
    }
    slab_iterator operator++(int) // postfix increment
    {
        auto copy = *this;
        ++(*this);
        return copy;
    }

    slab_iterator& operator--()
    {
        --slab_idx;
        return *this;
    }
    slab_iterator operator--(int) // postfix decrement
    {
        auto copy = *this;
        --(*this);
        return copy;
    }

    T span_begin() const
    {
        if constexpr (dim > 1)
            return box.lower[dim - 2];
        return 0;
    }
    T span_end() const // clamped, see span_iterator::_span_size
    {
        if constexpr (dim > 1)
            return std::max<T>(box.lower[dim - 2], box.upper[dim - 2] + 1);
        return 1;
    }

protected:
    Box<T, dim> const& box;
    T slab_idx;
};


/** @brief Iterator over the points of one span, compared against the span's end
 * coordinate along the fastest-varying dimension.
 */
template<typename T, std::size_t dim>
class BoxSpan<T, dim>::point_iterator
{
public:
    point_iterator(Point<T, dim> const& start)
        : point{start}
    {
    }

    auto& operator*() const { return point; }

    point_iterator& operator++()
    {
        ++point[dim - 1];
        return *this;
    }

    bool operator!=(T const end) const { return point[dim - 1] != end; }

private:
    Point<T, dim> point;
};


template<typename T, std::size_t dim>
auto make_box_slabs(Box<T, dim> const box)
{
    return BoxSlabs<T, dim>{box};
}


/** @brief Calls fn(point, args...) for every point of `box`, walking it span-wise.
 *
 * `point` is owned by the span's point_iterator: it is only valid for the duration
 * of the call.
 * `args` are passed to every call as lvalues, never forwarded: fn runs once per
 * point, so a by-value parameter could otherwise move from them on the first one.
 */
template<typename T, std::size_t dim>
void for_each_point(Box<T, dim> const& box, auto&& fn, auto&&... args)
{
    for (auto const& slab : make_box_slabs(box))
        for (auto const& span : slab)
            for (auto const& point : span)
                fn(point, args...);
}


} // namespace PHARE::core


#endif //  PHARE_CORE_UTILITIES_BOX_BOX_SLABS_HPP
