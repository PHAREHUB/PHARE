#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP


#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core
{

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxRows
{
    using raw_value_type = Array_t::value_type;
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    Array_t& arr;
    Box<std::uint32_t, dim> const& box;
    std::uint32_t const slab_idx;
    std::uint32_t row_idx;
    std::uint32_t const row_size = _row_size();
    Span<value_type> row{0, 0};

    FieldBoxRows& operator++()
    {
        ++row_idx;
        return *this;
    }

    auto& operator*()
    {
        if constexpr (dim == 1)
            return (row = Span<value_type>{&arr(box.lower[0]), row_size});
        if constexpr (dim == 2)
            return (row = Span<value_type>{&arr(row_idx, box.lower[1]), row_size});
        if constexpr (dim == 3)
            return (row = Span<value_type>{&arr(slab_idx, row_idx, box.lower[2]), row_size});
    }

    bool operator==(FieldBoxRows const& that) const { return row_idx == that.row_idx; }
    bool operator!=(FieldBoxRows const& that) const { return row_idx != that.row_idx; }

    auto _row_size() const { return box.upper[dim - 1] - box.lower[dim - 1] + 1; }
};


template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxSlab
{
    using FieldBoxRows_t = FieldBoxRows<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> const& box;
    std::uint32_t slab_idx;

    FieldBoxRows_t begin() { return {arr, box, slab_idx, row_begin()}; }
    FieldBoxRows_t begin() const { return {arr, box, slab_idx, row_begin()}; }
    FieldBoxRows_t end() { return {arr, box, slab_idx, row_end()}; }
    FieldBoxRows_t end() const { return {arr, box, slab_idx, row_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

    bool operator==(FieldBoxSlab const& that) const { return slab_idx == that.slab_idx; }
    bool operator!=(FieldBoxSlab const& that) const { return slab_idx != that.slab_idx; }

    FieldBoxSlab& operator++()
    {
        ++slab_idx;
        return *this;
    }

    std::uint32_t row_begin() const
    {
        if constexpr (dim > 1)
            return box.lower[dim - 2];
        return 0;
    }
    std::uint32_t row_end() const
    {
        if constexpr (dim > 1)
            return box.upper[dim - 2] + 1;
        return 1;
    }
};

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxSpan
{
    using FieldBoxSlab_t = FieldBoxSlab<Array_t>;

    Box<std::uint32_t, dim> const box;
    Array_t& arr;

    FieldBoxSlab_t begin() { return {arr, box, slab_begin()}; }
    FieldBoxSlab_t begin() const { return {arr, box, slab_begin()}; }
    FieldBoxSlab_t end() { return {arr, box, slab_end()}; }
    FieldBoxSlab_t end() const { return {arr, box, slab_end()}; }

    std::uint32_t slab_begin() const
    {
        if constexpr (dim > 2)
            return box.lower[0];
        return 0;
    }
    std::uint32_t slab_end() const
    {
        if constexpr (dim > 2)
            return box.upper[0] + 1;
        return 1;
    }
};


template<typename Array_t, std::size_t dim>
auto make_field_box_span(Box<std::uint32_t, dim> const box, Array_t& arr)
{
    return FieldBoxSpan<Array_t>{box, arr};
}

template<typename Array_t, std::size_t dim>
auto make_field_box_span(Box<std::uint32_t, dim> const box, Array_t const& arr)
{
    return FieldBoxSpan<Array_t const>{box, arr};
}

} // namespace PHARE::core


#endif //  PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
