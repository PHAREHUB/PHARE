// C++11 version of https://en.cppreference.com/w/cpp/container/span

#ifndef PHARE_CORE_UTILITIES_SPAN_H
#define PHARE_CORE_UTILITIES_SPAN_H

#include <cstddef> // for size_t
#include <numeric> // for accumulate
#include <vector>  // for vector
#include "core/utilities/types.h"

namespace PHARE::core
{
template<typename T, typename SIZE = size_t>
struct Span
{
    using value_type = T;

    T& operator[](SIZE i) const { return ptr[i]; }
    T const* const& data() const { return ptr; }
    T const* const& begin() const { return ptr; }
    T* end() const { return ptr + s; }
    SIZE const& size() const { return s; }


    T const* ptr = nullptr;
    SIZE s       = 0;
};

template<typename T, typename SIZE = size_t>
struct SpanSet
{
    using value_type = T;
    using SpanSet_   = SpanSet<T, SIZE>;

    SpanSet() = default;

    SpanSet(std::vector<SIZE>&& sizes_)
        : size{std::accumulate(sizes_.begin(), sizes_.end(), 0)}
        , sizes(sizes_)
        , displs(core::displacementFrom(sizes))
        , vec(size)
    {
    }

    SpanSet(SpanSet&& from)
        : size{from.size}
        , sizes(std::move(from.sizes))
        , displs(std::move(from.displs))
        , vec(std::move(from.vec))
    {
    }

    Span<T, SIZE> operator[](SIZE i) const
    {
        return {this->vec.data() + displs[i], this->sizes[i]};
    }

    T* data() const { return const_cast<T*>(vec.data()); }

    struct iterator
    {
        iterator(SpanSet_* _sv)
            : sv(_sv)
        {
        }
        iterator operator++()
        {
            curr_pos += sv->sizes[curr_ptr++];
            return *this;
        }
        bool operator!=(iterator const& other) const { return curr_ptr != sv->sizes.size(); }
        Span<T, SIZE> operator*() const { return {sv->vec.data() + curr_pos, sv->sizes[curr_ptr]}; }

        SpanSet_* sv  = nullptr;
        SIZE curr_pos = 0, curr_ptr = 0;
    };

    auto begin() { return iterator(this); }
    auto cbegin() const { return iterator(this); }

    auto end() { return iterator(this); }
    auto cend() const { return iterator(this); }

    SIZE size;
    std::vector<SIZE> sizes;
    std::vector<SIZE> displs;
    std::vector<T> vec;
};
} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_SPAN_H
