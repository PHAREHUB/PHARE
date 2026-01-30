// C++11 version of https://en.cppreference.com/w/cpp/container/span

#ifndef PHARE_CORE_UTILITIES_SPAN_HPP
#define PHARE_CORE_UTILITIES_SPAN_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"

#include <vector>
#include <cstddef>
#include <numeric>

namespace PHARE::core
{

template<typename T>
concept Spannable = requires(T t) {
    { t.size() };
    { t.data() };
};


template<typename T, typename SIZE = std::size_t>

struct Span
{
    using value_type = T;

    NO_DISCARD auto& operator[](SIZE i) { return ptr[i]; }
    NO_DISCARD auto& operator[](SIZE i) const { return ptr[i]; }
    NO_DISCARD T const* const& data() const { return ptr; }
    NO_DISCARD T const* const& begin() const { return ptr; }
    NO_DISCARD T* end() const { return ptr + s; }
    NO_DISCARD SIZE const& size() const { return s; }

    T const* ptr = nullptr;
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

    NO_DISCARD Span<T, SIZE> operator[](SIZE i) const
    {
        return {this->vec.data() + displs[i], this->sizes[i]};
    }

    NO_DISCARD T* data() const { return const_cast<T*>(vec.data()); }

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

    NO_DISCARD auto begin() { return iterator(this); }
    NO_DISCARD auto cbegin() const { return iterator(this); }

    NO_DISCARD auto end() { return iterator(this); }
    NO_DISCARD auto cend() const { return iterator(this); }

    SIZE size;
    std::vector<SIZE> sizes;
    std::vector<SIZE> displs;
    std::vector<T> vec;
};
} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_SPAN_HPP
