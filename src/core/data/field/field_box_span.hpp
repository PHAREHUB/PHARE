#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP


#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core
{

template<typename Array_t>
class FieldBoxSpans
{
    using raw_value_type = Array_t::value_type;

    constexpr auto static dim = Array_t::dimension;

public:
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    FieldBoxSpans(Array_t& arr, Box<std::uint32_t, dim> const& box, std::uint32_t const slab_idx,
                  std::uint32_t const span_idx)
        : arr{arr}
        , box{box}
        , slab_idx{slab_idx}
        , span_idx{span_idx}
        , span_size{box.upper[dim - 1] - box.lower[dim - 1] + 1}
    {
    }

    FieldBoxSpans& operator++()
    {
        ++span_idx;
        return *this;
    }

    auto& operator*()
    {
        if constexpr (dim == 1)
            return (span = Span<value_type>{&arr(box.lower[0]), span_size});
        if constexpr (dim == 2)
            return (span = Span<value_type>{&arr(span_idx, box.lower[1]), span_size});
        if constexpr (dim == 3)
            return (span = Span<value_type>{&arr(slab_idx, span_idx, box.lower[2]), span_size});
    }

    bool operator==(FieldBoxSpans const& that) const { return span_idx == that.span_idx; }
    bool operator!=(FieldBoxSpans const& that) const { return span_idx != that.span_idx; }


private:
    Array_t& arr;
    Box<std::uint32_t, dim> const& box;
    std::uint32_t const slab_idx;
    std::uint32_t span_idx;
    std::uint32_t const span_size;
    Span<value_type> span{0, 0};
};


template<typename Array_t>
class FieldBoxSlab
{
    using FieldBoxSpans_t     = FieldBoxSpans<Array_t>;
    constexpr auto static dim = Array_t::dimension;

public:
    FieldBoxSlab(Array_t& arr, Box<std::uint32_t, dim> const& box, std::uint32_t const slab_idx)
        : arr{arr}
        , box{box}
        , slab_idx{slab_idx}
    {
    }

    FieldBoxSpans_t begin() { return {arr, box, slab_idx, spans_begin()}; }
    FieldBoxSpans_t begin() const { return {arr, box, slab_idx, spans_begin()}; }
    FieldBoxSpans_t end() { return {arr, box, slab_idx, spans_end()}; }
    FieldBoxSpans_t end() const { return {arr, box, slab_idx, spans_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

    bool operator==(FieldBoxSlab const& that) const { return slab_idx == that.slab_idx; }
    bool operator!=(FieldBoxSlab const& that) const { return slab_idx != that.slab_idx; }

    FieldBoxSlab& operator++()
    {
        ++slab_idx;
        return *this;
    }

    std::uint32_t spans_begin() const
    {
        if constexpr (dim > 1)
            return box.lower[dim - 2];
        return 0;
    }
    std::uint32_t spans_end() const
    {
        if constexpr (dim > 1)
            return box.upper[dim - 2] + 1;
        return 1;
    }

private:
    Array_t& arr;
    Box<std::uint32_t, dim> const& box;
    std::uint32_t slab_idx;
};

template<typename Array_t>
class FieldBoxSpan
{
    using FieldBoxSlab_t      = FieldBoxSlab<Array_t>;
    constexpr auto static dim = Array_t::dimension;

public:
    FieldBoxSpan(Array_t& arr, Box<std::uint32_t, dim> const& box)
        : arr{arr}
        , box{box}
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
    Box<std::uint32_t, dim> const box;
};


template<typename Array_t, std::size_t dim>
auto make_field_box_span(Box<std::uint32_t, dim> const& box, Array_t& arr)
{
    return FieldBoxSpan<Array_t>{arr, box};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_span(Box<std::uint32_t, dim> const& box, Array_t const& arr)
{
    return FieldBoxSpan<Array_t const>{arr, box};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
