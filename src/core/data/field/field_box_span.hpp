#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP


#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/box/box_span.hpp"

#include <cstddef>
#include <cstdint>
#include <tuple>


namespace PHARE::core
{

template<typename Array_t>
struct FieldBoxRows : BoxRows<std::uint32_t, Array_t::dimension>
{
    auto constexpr static dim = Array_t::dimension;
    using Super               = BoxRows<std::uint32_t, dim>;
    using Super::point;
    using Super::row_size;
    using raw_value_type = Array_t::value_type;
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    FieldBoxRows(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    Array_t& arr;
    Span<value_type> row{0, 0};

    auto& operator*() { return (row = Span<value_type>{&arr(point()), row_size}); }
};


template<typename Array_t>
struct FieldBoxPointRows : public FieldBoxRows<Array_t>
{
    using Super               = FieldBoxRows<Array_t>;
    using value_type          = Super::value_type;
    auto constexpr static dim = Array_t::dimension;

    FieldBoxPointRows(auto&&... args)
        : Super{args...}
        , tup{std::forward_as_tuple(*super(), Super::_point)}
    {
    }

    Super& super() { return *this; }
    auto& operator*()
    {
        std::get<0>(tup) = *super();
        std::get<1>(tup) = Super::point();
        return tup;
    }

    std::tuple<Span<value_type>&, Point<std::uint32_t, dim>&> tup;
};




template<typename Array_t, typename Rows_t = FieldBoxRows<Array_t>>
struct FieldBoxSlab : BoxSlab<std::uint32_t, Array_t::dimension>
{
    auto constexpr static dim = Array_t::dimension;
    using FieldBoxRows_t      = Rows_t;
    using BaseBoxRows_t       = BoxRows<std::uint32_t, Array_t::dimension>;
    using Super               = BoxSlab<std::uint32_t, dim>;
    using Super::box;
    using Super::row_begin;
    using Super::row_end;
    using Super::slab_idx;


    FieldBoxSlab(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    Array_t& arr;

    FieldBoxRows_t begin() { return {arr, box, slab_idx, row_begin()}; }
    FieldBoxRows_t begin() const { return {arr, box, slab_idx, row_begin()}; }
    BaseBoxRows_t end() { return {box, slab_idx, row_end()}; }
    BaseBoxRows_t end() const { return {box, slab_idx, row_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }
};

template<typename Array_t, typename Rows_t = FieldBoxRows<Array_t>>
struct FieldBoxSpan : public BoxSpan<std::uint32_t, Array_t::dimension>
{
    auto constexpr static dim = Array_t::dimension;
    using FieldBoxSlab_t      = FieldBoxSlab<Array_t, Rows_t>;
    using Super               = BoxSpan<std::uint32_t, dim>;
    using Super::box;
    using Super::slab_begin;
    using Super::slab_end;

    FieldBoxSpan(Array_t& ar, auto&&... args)
        : Super{args...}
        , arr{ar}
    {
    }

    Array_t& arr;

    FieldBoxSlab_t begin() { return {arr, box, slab_begin()}; }
    FieldBoxSlab_t begin() const { return {arr, box, slab_begin()}; }
    FieldBoxSlab_t end() { return {arr, box, slab_end()}; }
    FieldBoxSlab_t end() const { return {arr, box, slab_end()}; }
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
    using Row_t = FieldBoxPointRows<Array_t>;
    return FieldBoxSpan<Array_t, Row_t>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_point_span(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    using Row_t = FieldBoxPointRows<Array_t const>;
    return FieldBoxSpan<Array_t const, Row_t>{arr, box};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
