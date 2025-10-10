#ifndef CORE_NUMERICS_SLOPE_LIMITER_MIN_MOD_HPP
#define CORE_NUMERICS_SLOPE_LIMITER_MIN_MOD_HPP

#include <algorithm>
#include <cstdlib>
#include <type_traits>

namespace PHARE::core
{
struct MinModLimiter
{
    template<typename T, typename... Args>
    static auto limit(T const& first, Args const&... rest)
    {
        bool all_positive = (first > 0) && ((rest > 0) && ...);
        bool all_negative = (first < 0) && ((rest < 0) && ...);

        if (!all_positive && !all_negative)
        {
            return static_cast<std::common_type_t<T, Args...>>(0);
        }

        auto min_abs = [](auto a, auto b) { return std::abs(a) < std::abs(b); };

        return std::min({first, rest...}, min_abs);
    }
};
} // namespace PHARE::core

#endif
