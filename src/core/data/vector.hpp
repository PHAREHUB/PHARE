#ifndef PHARE_CORE_DATA_VECTOR_HPP
#define PHARE_CORE_DATA_VECTOR_HPP



#include <vector>
#include <cstdint>

namespace PHARE::core
{

// Will automagically shrink itself if the requested sizes are
//  below a threshold for some number of requests
template<typename T>
struct MinimizingVector
{
    using vector_t = std::vector<T>;

    // resize to s, preserving existing content
    auto& get(std::size_t s)
    {
        _maybe_shrink<true>(s);
        vec.resize(s);
        return *this;
    }

    // resize to s, discarding existing content
    auto& get_no_copy(std::size_t s)
    {
        _maybe_shrink<false>(s);
        vec.resize(s);
        return *this;
    }

    // reserve capacity for s elements, size stays 0
    auto& reserve_and_clear(std::size_t s)
    {
        _maybe_shrink<false>(s);
        vec.clear();
        vec.reserve(s);
        return *this;
    }

    auto& resize(std::size_t s) { return get(s); }
    auto& clear()
    {
        vec.clear();
        return *this;
    }
    void destroy() { vec = vector_t{}; }

    auto& operator[](std::size_t i) { return vec.data()[i]; }
    auto& operator[](std::size_t i) const { return vec.data()[i]; }
    auto& operator*() { return vec; }
    auto& operator()() { return vec; }
    auto& operator()() const { return vec; }

    auto& emplace_back(auto&&... args)
    {
        return vec.emplace_back(std::forward<decltype(args)>(args)...);
    }
    void push_back(auto&&... args) { vec.push_back(std::forward<decltype(args)>(args)...); }
    auto data() { return vec.data(); }
    auto data() const { return vec.data(); }
    auto size() const { return vec.size(); }
    auto capacity() const { return vec.capacity(); }


    // public overridable
    double const percentile  = .80;
    double const realloc_to  = .85;
    std::size_t const period = 100;

private:
    template<bool copy_old>
    void _maybe_shrink(std::size_t s)
    {
        if (s < vec.capacity() * percentile)
            ++curr_access;
        else
            curr_access = 0;

        if (curr_access != period)
            return;

        auto const new_cap = static_cast<std::size_t>(vec.capacity() * realloc_to);
        if constexpr (copy_old)
        {
            vector_t r;
            r.reserve(new_cap);
            r   = vec;
            vec = std::move(r);
        }
        else
        {
            vec = vector_t{};
            vec.reserve(new_cap);
        }
        curr_access = 0;
    }

    vector_t vec{};
    std::uint16_t curr_access = 0;
};


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_VECTOR_HPP */
