#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SLABS_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SLABS_HPP


#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_slabs.hpp"

#include <cstddef>


namespace PHARE::core
{

/** @brief Like BoxSpans, but dereferencing yields a Span<value_type> over the
 * underlying array's data for the current span, instead of a (point, size) pair.
 */
template<typename Array_t>
class FieldBoxSpans : public BoxSpans<std::uint32_t, Array_t::dimension>
{
    constexpr auto static dim = Array_t::dimension;
    using Super               = BoxSpans<std::uint32_t, dim>;
    using raw_value_type      = Array_t::value_type;
    using Super::span_size;

protected:
    using Super::point;

public:
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    FieldBoxSpans(Array_t& arr, auto&&... args)
        : Super{args...}
        , arr{arr}
    {
    }


    auto& operator*() { return (span = Span<value_type>{&arr(point()), span_size}); }

private:
    Array_t& arr;
    Span<value_type> span{0, 0};
};


/** @brief Like FieldBoxSpans, but dereferencing also yields the current span's
 * starting point alongside its data Span, for callers that need to track the
 * box index in step with the field data (e.g. cross-level refinement, which
 * must know where each span sits on both the coarse and fine grids).
 */
template<typename Array_t>
class FieldBoxPointSpans : public FieldBoxSpans<Array_t>
{
    using Super = FieldBoxSpans<Array_t>;

public:
    using value_type = Super::value_type;
    using Super::Super;

    // a tuple of references to this iterator's own span and point, rebuilt per dereference so
    //  copies never alias the source's members
    auto operator*() { return std::forward_as_tuple(Super::operator*(), this->point()); }
};



template<typename Array_t, typename Spans_t = FieldBoxSpans<Array_t>>
class FieldBoxSlab : public BoxSlab<std::uint32_t, Array_t::dimension>
{
    auto constexpr static dim = Array_t::dimension;
    using FieldBoxSpans_t     = Spans_t;
    using BaseBoxSpans_t      = BoxSpans<std::uint32_t, Array_t::dimension>;
    using Super               = BoxSlab<std::uint32_t, dim>;
    using Super::box;
    using Super::slab_idx;
    using Super::span_begin;
    using Super::span_end;


public:
    FieldBoxSlab(Array_t& ar, auto&&... args)
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



template<typename Array_t, typename Spans_t = FieldBoxSpans<Array_t>>
class FieldBoxSlabs : public BoxSlabs<std::uint32_t, Array_t::dimension>
{
    constexpr auto static dim = Array_t::dimension;
    using FieldBoxSlab_t      = FieldBoxSlab<Array_t, Spans_t>;
    using Super               = BoxSlabs<std::uint32_t, dim>;
    using Super::box;
    using Super::slab_begin;
    using Super::slab_end;

public:
    FieldBoxSlabs(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    FieldBoxSlab_t begin() { return {arr, box, slab_begin()}; }
    FieldBoxSlab_t begin() const { return {arr, box, slab_begin()}; }
    FieldBoxSlab_t end() { return {arr, box, slab_end()}; }
    FieldBoxSlab_t end() const { return {arr, box, slab_end()}; }

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
    using PointSpans_t = FieldBoxPointSpans<Array_t>;
    return FieldBoxSlabs<Array_t, PointSpans_t>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_point_slabs(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    using PointSpans_t = FieldBoxPointSpans<Array_t const>;
    return FieldBoxSlabs<Array_t const, PointSpans_t>{arr, box};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_DATA_FIELD_FIELD_BOX_SLABS_HPP
