#ifndef PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/box/box.hpp"

#include <cstddef>


namespace PHARE::core // 3d only for now
{

template<std::size_t dim>
struct BoxRow
{
    Point<std::uint32_t, dim> start_;
    std::uint32_t size_;

    auto& local_start() const { return start_; }
    auto& size() const { return size_; }
};




template<std::size_t dim>
struct BoxSpan
{
    std::uint32_t constexpr static ZERO = 0;

    template<typename Span_t>
    struct BoxRows
    {
        Span_t& span;

        std::uint32_t k;
        std::uint32_t j;
        std::uint32_t s = span._size_i();
        BoxRow<dim> iter{begin(), s};

        Point<std::uint32_t, dim> begin() const
        {
            if constexpr (dim == 1)
                return Point<std::uint32_t, dim>{ZERO};

            if constexpr (dim == 2)
                return Point<std::uint32_t, dim>{j, ZERO};

            if constexpr (dim == 3)
                return Point<std::uint32_t, dim>{k, j, ZERO};
        }

        BoxRows& operator++()
        {
            ++j;
            if constexpr (dim == 1)
                iter.start_ = {span.box.lower[0]};
            if constexpr (dim == 2)
                iter.start_ = {j, span.box.lower[1]};
            if constexpr (dim == 3)
                iter.start_ = {k, j, span.box.lower[2]};
            // PHARE_LOG_LINE_SS(k);
            return *this;
        }

        auto& operator*() const { return iter; }
        bool operator==(BoxRows const& that) const { return j == that.j; }
        bool operator!=(BoxRows const& that) const { return j != that.j; }

        void next()
        {
            j = 0;
            ++k;
            if constexpr (dim == 1)
                iter.start_ = Point{span.box.lower[0]};
            if constexpr (dim == 2)
                iter.start_ = Point{j, span.box.lower[1]};
            if constexpr (dim == 3)
                iter.start_ = Point{k, j, span.box.lower[2]};
            // PHARE_LOG_LINE_SS(k);
        }
    };

    template<typename Span_t>
    struct BoxSlab
    {
        Span_t& span;
        BoxRows<Span_t> br{span, 0, 0};

        BoxRows<Span_t> begin() const { return br; }
        BoxRows<Span_t> end() const { return {span, span._end_k(), span._end_j()}; }

        void next() { br.next(); }
    };

    template<typename Span_t>
    struct BoxSlabber
    {
        Span_t& span;
        std::uint32_t k = 0;
        BoxSlab<Span_t> slab{span};

        bool operator==(BoxSlabber const& that) const { return k == that.k; }
        bool operator!=(BoxSlabber const& that) const { return k != that.k; }

        BoxSlabber& operator++()
        {
            slab.next();
            ++k;
            return *this;
        }
        auto& operator*() { return slab; }
    };

    Box<std::uint32_t, dim> const box;

    auto begin() { return BoxSlabber<BoxSpan<dim>>{*this}; }
    auto end() { return BoxSlabber<BoxSpan<dim>>{*this, _end_k()}; }
    auto begin() const { return BoxSlabber<BoxSpan<dim> const>{*this}; }
    auto end() const { return BoxSlabber<BoxSpan<dim> const>{*this, _end_k()}; }

    std::uint32_t _end_k() const
    {
        if constexpr (dim == 3)
            return box.upper[0] + 1;
        else
            return 1;
    }

    std::uint32_t _end_j() const
    {
        if constexpr (dim > 1)
            return box.upper[1] + 1;
        else
            return 1;
    }

    std::uint32_t _size_i() const { return box.upper[dim - 1] + 1; }
};



template<std::size_t dim>
auto make_box_span(Box<std::uint32_t, dim> const box)
{
    return BoxSpan<dim>{box};
}


} // namespace PHARE::core


#endif //  PHARE_CORE_UTILITIES_BOX_BOX_SPAN_HPP
