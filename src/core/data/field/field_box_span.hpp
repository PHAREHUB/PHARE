#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP


#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_span.hpp"

#include <cstddef>


namespace PHARE::core
{

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


template<typename Array_t>
class FieldBoxPointSpans : public FieldBoxSpans<Array_t>
{
    auto constexpr static dim = Array_t::dimension;
    using Super               = FieldBoxSpans<Array_t>;

public:
    using value_type = Super::value_type;

    FieldBoxPointSpans(auto&&... args)
        : Super{args...}
        , tup{std::forward_as_tuple(*super(), Super::_point)}
    {
    }

    Super& super() { return *this; }
    auto& operator*()
    {
        std::get<0>(tup) = *super();
        std::get<1>(tup) = this->point();
        return tup;
    }

private:
    std::tuple<Span<value_type>&, Point<std::uint32_t, dim>&> tup;
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
class FieldBoxSpan : public BoxSpan<std::uint32_t, Array_t::dimension>
{
    constexpr auto static dim = Array_t::dimension;
    using FieldBoxSlab_t      = FieldBoxSlab<Array_t, Spans_t>;
    using Super               = BoxSpan<std::uint32_t, dim>;
    using Super::box;
    using Super::slab_begin;
    using Super::slab_end;

public:
    FieldBoxSpan(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    FieldBoxSlab_t begin() { return {arr, box, slabs_begin()}; }
    FieldBoxSlab_t begin() const { return {arr, box, slabs_begin()}; }
    FieldBoxSlab_t end() { return {arr, box, slabs_end()}; }
    FieldBoxSlab_t end() const { return {arr, box, slabs_end()}; }

    std::uint32_t slabs_begin() const
    {
        if constexpr (dim > 2)
            return box.lower[0];
        return 0;
    }
    std::uint32_t slabs_end() const
    {
        if constexpr (dim > 2)
            return box.upper[0] + 1;
        return 1;
    }

private:
    Array_t& arr;
};



template<typename Array_t, std::size_t dim>
auto make_field_box_span(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    return FieldBoxSpan<Array_t>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_span(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    return FieldBoxSpan<Array_t const>{arr, box};
}


template<typename Array_t, std::size_t dim>
auto make_field_box_point_span(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    using Row_t = FieldBoxPointSpans<Array_t>;
    return FieldBoxSpan<Array_t, Row_t>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_point_span(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    using Row_t = FieldBoxPointSpans<Array_t const>;
    return FieldBoxSpan<Array_t const, Row_t>{arr, box};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
