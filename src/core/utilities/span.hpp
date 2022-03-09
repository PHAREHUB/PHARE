// C++11 version of https://en.cppreference.com/w/cpp/container/span

#ifndef PHARE_CORE_UTILITIES_SPAN_HPP
#define PHARE_CORE_UTILITIES_SPAN_HPP

#include <vector>
#include <cstddef>
#include <numeric>

#include "core/logger.hpp"
#include "core/utilities/types.hpp"

namespace PHARE::core
{
template<typename T, typename SIZE = std::size_t>
struct Span
{
    using value_type = std::decay_t<T>;

    auto& operator[](SIZE i) { return ptr[i]; }
    auto& operator[](SIZE i) const { return ptr[i]; }
    auto data() const { return ptr; }
    auto data() { return ptr; }
    auto begin() const { return ptr; }
    auto end() const { return ptr + s; }
    SIZE const& size() const { return s; }

    T* const ptr = nullptr;
    SIZE s       = 0;
};


template<typename T, typename SIZE = std::size_t>
class VectorSpan : private StackVar<std::vector<T>>, public core::Span<T, SIZE>
{
    using Vector = StackVar<std::vector<T>>;
    using Span_  = Span<T, SIZE>;

public:
    VectorSpan(std::size_t size, T value)
        : Vector{std::vector<T>(size, value)}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
    VectorSpan(std::vector<T>&& vec_)
        : Vector{std::move(vec_)}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
    VectorSpan(std::vector<T> const& vec_)
        : Vector{vec_}
        , Span_{Vector::var.data(), Vector::var.size()}
    {
    }
};



template<typename T, typename SIZE = std::size_t>
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

    Span<T, SIZE> operator[](SIZE i) { return {this->vec.data() + displs[i], this->sizes[i]}; }
    Span<T, SIZE> operator[](SIZE i) const
    {
        return {this->vec.data() + displs[i], this->sizes[i]};
    }

    T* data() const { return const_cast<T*>(vec.data()); }
    T* data() { return const_cast<T*>(vec.data()); }

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




template<typename T, std::size_t size>
auto flatten(std::vector<std::array<T, size>> const& data)
{
    assert(data.size() > 0);

    return Span<T const, std::size_t>{data.data()->data(), data.size() * size};
}


template<typename T, std::size_t size>
auto flatten(std::vector<std::array<T, size>>& data)
{
    assert(data.size() > 0);

    return Span<T, std::size_t>{data.data()->data(), data.size() * size};
}



} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_SPAN_HPP
