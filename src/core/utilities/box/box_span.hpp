#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP


// #include "core/def.hpp"
// #include "core/logger.hpp"
// #include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core
{

template<typename T, std::size_t dim_>
class BoxSpans
{
    auto constexpr static dim = dim_;

public:
    using value_type = T;

    BoxSpans(Box<T, dim> const& box, T const slab_idx, T const span_idx)
        : box{box}
        , slab_idx{slab_idx}
        , span_idx{span_idx}
        , span_size{_span_size()}
    {
    }

    BoxSpans& operator++()
    {
        ++span_idx;
        return *this;
    }

    auto& point()
    {
        if constexpr (dim == 1)
            return (_point = {box.lower[0]});
        if constexpr (dim == 2)
            return (_point = {span_idx, box.lower[1]});
        if constexpr (dim == 3)
            return (_point = {slab_idx, span_idx, box.lower[2]});
    }

    auto operator*() { return std::forward_as_tuple(point(), span_size); }

    bool operator==(BoxSpans const& that) const { return span_idx == that.span_idx; }
    bool operator!=(BoxSpans const& that) const { return span_idx != that.span_idx; }


protected:
    std::uint32_t _span_size() const { return box.upper[dim - 1] - box.lower[dim - 1] + 1; }

    Box<T, dim> const& box;
    T const slab_idx;
    T span_idx;
    std::uint32_t const span_size;
    Point<T, dim> _point{};
};


template<typename T, std::size_t dim_>
class BoxSlab
{
    auto constexpr static dim = dim_;
    using BoxSpans_t          = BoxSpans<T, dim>;

public:
    using value_type = T;

    BoxSlab(Box<T, dim> const& _box, T const _slab_idx)
        : box{_box}
        , slab_idx{_slab_idx}
    {
    }

    BoxSpans_t begin() { return {box, slab_idx, span_begin()}; }
    BoxSpans_t begin() const { return {box, slab_idx, span_begin()}; }
    BoxSpans_t end() { return {box, slab_idx, span_end()}; }
    BoxSpans_t end() const { return {box, slab_idx, span_end()}; }

    auto& operator*() { return *this; }
    auto& operator*() const { return *this; }

    bool operator==(BoxSlab const& that) const { return slab_idx == that.slab_idx; }
    bool operator!=(BoxSlab const& that) const { return slab_idx != that.slab_idx; }

    BoxSlab& operator++()
    {
        ++slab_idx;
        return *this;
    }

    T span_begin() const
    {
        if constexpr (dim > 1)
            return box.lower[dim - 2];
        return 0;
    }
    T span_end() const
    {
        if constexpr (dim > 1)
            return box.upper[dim - 2] + 1;
        return 1;
    }

protected:
    Box<T, dim> const& box;
    T slab_idx;
};

template<typename T, std::size_t dim_>
class BoxSpan
{
    auto constexpr static dim = dim_;
    using BoxSlab_t           = BoxSlab<T, dim>;

public:
    using value_type = T;

    BoxSpan(Box<T, dim> const& _box)
        : box{_box}
    {
    }

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

protected:
    Box<T, dim> const box;
};


template<typename T, std::size_t dim>
auto make_box_span(Box<T, dim> const box)
{
    return BoxSpan<T, dim>{box};
}


} // namespace PHARE::core


#endif //  PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
