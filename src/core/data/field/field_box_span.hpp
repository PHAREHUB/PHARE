#ifndef PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_BOX_SPAN_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core // 3d only for now
{

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxRows
{
    using raw_value_type = Array_t::value_type;
    using value_type
        = std::conditional_t<std::is_const_v<Array_t>, raw_value_type const, raw_value_type>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    std::uint32_t k;
    std::uint32_t j = box.lower[1];
    std::uint32_t s = box.upper[2] - box.lower[2] + 1;


    Span<value_type> row{&arr(k, j, box.lower[2]), s};

    FieldBoxRows& operator++()
    {
        row = Span<value_type>{&arr(k, ++j, box.lower[2]), s};
        return *this;
    }
    auto& operator*() { return row; }

    bool operator==(FieldBoxRows const& that) const { return j == that.j; }
    bool operator!=(FieldBoxRows const& that) const { return j != that.j; }

    auto point() const { return Point{k, j, box.lower[2]}; }

    void next()
    {
        j   = box.lower[1];
        row = Span<value_type>{&arr(++k, j, box.lower[2]), s};
    }
};

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxSlab
{
    using FieldBoxRows_t = FieldBoxRows<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    bool _end = false;
    FieldBoxRows_t br{arr, box, _end ? box.upper[0] : box.lower[0]};

    FieldBoxRows_t begin() const
    {
        assert(!_end);
        return br;
    }
    FieldBoxRows_t end() const { return {arr, box, box.upper[0], box.upper[1] /*+ 1*/, 0, {0, 0}}; }

    void next() { br.next(); }
};

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxSlabber
{
    using FieldBoxSlab_t = FieldBoxSlab<Array_t>;

    Array_t& arr;
    Box<std::uint32_t, dim> box;
    bool _end       = false;
    std::uint32_t k = _end ? box.upper[0] /*+ 1 */ : box.lower[0];

    FieldBoxSlab_t slab{arr, box, _end};

    bool operator==(FieldBoxSlabber const& that) const { return k == that.k; }
    bool operator!=(FieldBoxSlabber const& that) const { return k != that.k; }

    FieldBoxSlabber& operator++()
    {
        slab.next();
        ++k;
        return *this;
    }
    auto& operator*() { return slab; }
};

template<typename Array_t, std::size_t dim = Array_t::dimension>
struct FieldBoxSpan
{
    using FieldBoxSlabber_t = FieldBoxSlabber<Array_t>;

    Box<std::uint32_t, dim> box;
    Array_t& arr;

    FieldBoxSlabber_t b{arr, box};
    FieldBoxSlabber_t e{arr, box, true};

    auto begin() { return b; }
    auto end() { return e; }
    auto begin() const { return b; }
    auto end() const { return e; }
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
