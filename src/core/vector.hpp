#ifndef PHARE_CORE_VECTOR_HPP
#define PHARE_CORE_VECTOR_HPP

#include <vector>
#include <cstdint>
#include "core/utilities/span.hpp"

namespace PHARE::core
{


template<typename T, bool copy_old_ = true>
struct MinimizingVector
{
    template<bool copy_old = copy_old_>
    auto& get(std::size_t const& s)
    {
        if (s < v.capacity() * percentile)
            ++_c;
        else
            _c = 0;

        if (_c == period)
        {
            std::vector<T> r(v.capacity() * realloc_to);
            if constexpr (copy_old)
                r = v;
            v  = std::move(r);
            _c = 0;
        }

        v.resize(s);
        return v;
    }

    auto size() const { return v.size(); }
    auto capacity() const { return v.capacity(); }
    auto& operator()() { return v; }
    auto& operator()() const { return v; }
    auto& operator()(std::size_t const& s) { return get(s); }

    double const percentile  = .80;
    double const realloc_to  = .90;
    std::size_t const period = 100;

    std::vector<T> v{};
    std::uint16_t _c = 0;
};




} // namespace PHARE::core

#endif // PHARE_CORE_VECTOR_HPP
