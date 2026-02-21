#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core
{

template<typename T, std::size_t dim_>
struct BoxRows
{
    auto constexpr static dim = dim_;

    Box<T, dim> const& box;
    T const slab_idx;
    T row_idx;
    std::uint32_t const row_size = _row_size();
    Point<T, dim> _point{};

    BoxRows& operator++()
    {
        ++row_idx;
        return *this;
    }

    auto& point()
    {
        if constexpr (dim == 1)
            return (_point = {box.lower[0]});
        if constexpr (dim == 2)
            return (_point = {row_idx, box.lower[1]});
        if constexpr (dim == 3)
            return (_point = {slab_idx, row_idx, box.lower[2]});
    }

    auto operator*() { return std::forward_as_tuple(point(), row_size); }

    bool operator==(BoxRows const& that) const { return row_idx == that.row_idx; }
    bool operator!=(BoxRows const& that) const { return row_idx != that.row_idx; }

    std::uint32_t _row_size() const { return box.upper[dim - 1] - box.lower[dim - 1] + 1; }
};


template<typename T, std::size_t dim_>
struct BoxSlab
{
    auto constexpr static dim = dim_;
    using BoxRows_t           = BoxRows<T, dim>;

    Box<T, dim> const& box;
    T slab_idx;

    BoxRows_t begin() { return {box, slab_idx, row_begin()}; }
    BoxRows_t begin() const { return {box, slab_idx, row_begin()}; }
    BoxRows_t end() { return {box, slab_idx, row_end()}; }
    BoxRows_t end() const { return {box, slab_idx, row_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

    bool operator==(BoxSlab const& that) const { return slab_idx == that.slab_idx; }
    bool operator!=(BoxSlab const& that) const { return slab_idx != that.slab_idx; }

    BoxSlab& operator++()
    {
        ++slab_idx;
        return *this;
    }

    T row_begin() const
    {
        if constexpr (dim > 1)
            return box.lower[dim - 2];
        return 0;
    }
    T row_end() const
    {
        if constexpr (dim > 1)
            return box.upper[dim - 2] + 1;
        return 1;
    }
};

template<typename T, std::size_t dim_>
struct BoxSpan
{
    auto constexpr static dim = dim_;
    using BoxSlab_t           = BoxSlab<T, dim>;

    Box<T, dim> const box;

    BoxSlab_t begin() { return {box, slab_begin()}; }
    BoxSlab_t begin() const { return {box, slab_begin()}; }
    BoxSlab_t end() { return {box, slab_end()}; }
    BoxSlab_t end() const { return {box, slab_end()}; }

    T slab_begin() const
    {
        if constexpr (dim > 2)
            return box.lower[0];
        return 0;
    }
    T slab_end() const
    {
        if constexpr (dim > 2)
            return box.upper[0] + 1;
        return 1;
    }

    std::uint32_t size() const
    {
        auto s = slab_end() - slab_begin();
        assert(s > 0);
        return s;
    }
};


template<typename T, std::size_t dim>
auto make_box_span(Box<T, dim> const box)
{
    return BoxSpan<T, dim>{box};
}


} // namespace PHARE::core


#endif //  PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
