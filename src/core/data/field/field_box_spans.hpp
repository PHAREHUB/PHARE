#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SPANS_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SPANS_HPP

#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_slabs.hpp"

#include <cstddef>

namespace PHARE::core
{

/** @brief Iterates the slabs of a field box like BoxSlabs, but its spans also see the
 * field's data: span_iterator yields a Span<value_type> over the current span,
 * point_span_iterator (with_point) additionally yields the span's starting point.
 * See make_field_box_slabs() and make_field_box_point_slabs().
 */
template<typename Array_t, bool with_point = false>
class FieldBoxSlabs : public BoxSlabs<std::uint32_t, Array_t::dimension>
{
    constexpr auto static dim = Array_t::dimension;
    using Super               = BoxSlabs<std::uint32_t, dim>;
    using Super::box;
    using Super::slab_begin;
    using Super::slab_end;

    // spans run along the last index, contiguous in memory only for C ordering
    static_assert(Array_t::is_c_ordered, "FieldBoxSlabs requires a C ordered array");

public:
    class slab_iterator;
    class span_iterator;
    class point_span_iterator;

    FieldBoxSlabs(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    slab_iterator begin() { return {arr, box, slab_begin()}; }
    slab_iterator begin() const { return {arr, box, slab_begin()}; }
    slab_iterator end() { return {arr, box, slab_end()}; }
    slab_iterator end() const { return {arr, box, slab_end()}; }

private:
    Array_t& arr;
};


/** @brief Like BoxSlabs::span_iterator, but dereferencing yields a Span<value_type> over the
 * underlying array's data for the current span, instead of a (start point, span length) tuple.
 */
template<typename Array_t, bool with_point>
class FieldBoxSlabs<Array_t, with_point>::span_iterator
    : public BoxSlabs<std::uint32_t, dim>::span_iterator
{
    using Super          = BoxSlabs<std::uint32_t, dim>::span_iterator;
    using raw_value_type = Array_t::value_type;
    using Super::span_size;

protected:
    using Super::point;

public:
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    span_iterator(Array_t& arr, auto&&... args)
        : Super{args...}
        , arr{arr}
    {
    }


    auto& operator*() { return (span = Span<value_type>{&arr(point()), span_size}); }

private:
    Array_t& arr;
    Span<value_type> span{0, 0};
};


/** @brief Like span_iterator, but dereferencing also yields the current span's
 * starting point alongside its data Span, for callers that need to track the
 * box index in step with the field data (e.g. cross-level refinement, which
 * must know where each span sits on both the coarse and fine grids).
 */
template<typename Array_t, bool with_point>
class FieldBoxSlabs<Array_t, with_point>::point_span_iterator : public span_iterator
{
    using Super = span_iterator;

public:
    using value_type = Super::value_type;
    using Super::Super;

    // a tuple of references to this iterator's own span and point, rebuilt per dereference so
    //  copies never alias the source's members
    auto operator*() { return std::forward_as_tuple(Super::operator*(), this->point()); }
};


template<typename Array_t, bool with_point>
class FieldBoxSlabs<Array_t, with_point>::slab_iterator
    : public BoxSlabs<std::uint32_t, dim>::slab_iterator
{
    using FieldBoxSpans_t = std::conditional_t<with_point, point_span_iterator, span_iterator>;
    using BaseBoxSpans_t  = BoxSlabs<std::uint32_t, dim>::span_iterator;
    using Super           = BoxSlabs<std::uint32_t, dim>::slab_iterator;
    using Super::box;
    using Super::slab_idx;
    using Super::span_begin;
    using Super::span_end;

public:
    slab_iterator(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    FieldBoxSpans_t begin() { return {arr, box, slab_idx, span_begin()}; }
    FieldBoxSpans_t begin() const { return {arr, box, slab_idx, span_begin()}; }
    BaseBoxSpans_t end() { return {box, slab_idx, span_end()}; }
    BaseBoxSpans_t end() const { return {box, slab_idx, span_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

private:
    Array_t& arr;
};



template<typename Array_t, std::size_t dim>
auto make_field_box_slabs(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    return FieldBoxSlabs<Array_t>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_slabs(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    return FieldBoxSlabs<Array_t const>{arr, box};
}


template<typename Array_t, std::size_t dim>
auto make_field_box_point_slabs(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    return FieldBoxSlabs<Array_t, /*with_point*/ true>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_point_slabs(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    return FieldBoxSlabs<Array_t const, /*with_point*/ true>{arr, box};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_DATA_FIELD_FIELD_BOX_SPANS_HPP
