#ifndef PHARE_CORE_UTILITIES_NUMERIC_HPP
#define PHARE_CORE_UTILITIES_NUMERIC_HPP

#include "core/def.hpp"

#include <array>
#include <cmath>
#include <numeric>
#include <functional>
#include <type_traits>

namespace PHARE::core
{
template<typename Container, typename Multiplies = typename Container::value_type>
NO_DISCARD Multiplies product(Container const& container, Multiplies mul = 1)
{
    return std::accumulate(container.begin(), container.end(), mul, std::multiplies<Multiplies>());
}

template<typename Container, typename Return = typename Container::value_type>
NO_DISCARD Return sum(Container const& container, Return r = 0)
{
    return std::accumulate(container.begin(), container.end(), r);
}

template<typename T>
NO_DISCARD T norm(std::array<T, 3> vec)
{
    auto squarreSum = std::inner_product(std::begin(vec), std::end(vec), std::begin(vec), 0.);
    return std::sqrt(squarreSum);
}


template<typename Container, typename F>
NO_DISCARD auto sum_from(Container&& container, F fn)
{
    using value_type  = std::decay_t<Container>::value_type;
    using return_type = std::decay_t<std::invoke_result_t<F const&, value_type const&>>;
    return_type sum   = 0;
    for (auto const& el : container)
        sum += fn(el);
    return sum;
}


template<typename SignedInt, typename UnsignedInt>
bool diff_sign_int_equals(SignedInt const& i0, UnsignedInt const& i1)
{
    static_assert(std::is_unsigned_v<UnsignedInt>);
    static_assert(std::is_signed_v<SignedInt>);
    static_assert(sizeof(UnsignedInt) >= sizeof(SignedInt), "Bad int comparison!");
    if (i0 < 0)
        return false;
    return static_cast<UnsignedInt>(i0) == i1;
}


template<typename Int0, typename Int1>
bool int_equals(Int0 const& i0, Int1 const& i1)
{
    if constexpr (std::is_same_v<Int0, Int1>)
        return i0 == i1;
    else
    {
        if constexpr (std::is_unsigned_v<Int0> and std::is_signed_v<Int1>)
            return diff_sign_int_equals(i1, i0);
        if constexpr (std::is_unsigned_v<Int1> and std::is_signed_v<Int0>)
            return diff_sign_int_equals(i0, i1);
    }
    // reaching here == compiler error
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_NUMERIC_HPP
