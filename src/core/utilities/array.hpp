#ifndef PHARE_CORE_UTILITIES_ARRAY_HPP
#define PHARE_CORE_UTILITIES_ARRAY_HPP

#include "core/def.hpp"
#include "core/utilities/meta/for_N.hpp"

#include <array>
#include <limits>
#include <vector>
#include <cassert>
#include <cstdint>
#include <utility>
#include <type_traits>

namespace PHARE::core
{
template<typename Type, std::size_t Size> // std::array::fill is only constexpr in C++20 ffs
constexpr void fill(Type value, std::array<Type, Size>& array)
{
    for (std::size_t i = 0; i < Size; i++)
        array[i] = value;
}

template<std::size_t To_Size, std::size_t From_Size, typename Type>
NO_DISCARD constexpr std::array<Type, To_Size> sized_array(std::array<Type, From_Size> const& from)
{
    static_assert(To_Size <= From_Size, "invalid sized_array Size template, too large");

    if constexpr (To_Size == From_Size)
        return from;

    std::array<Type, To_Size> to{};

    for (std::size_t i = 0; i < to.size(); i++)
        to[i] = from[i];

    return to;
}


template<std::size_t To_Size, typename... Args>
NO_DISCARD constexpr auto as_sized_array(Args&&... args)
{
    auto arr = std::array{std::forward<decltype(args)>(args)...};

    return sized_array<To_Size>(arr);
}

template<typename Type, std::size_t size>
NO_DISCARD constexpr std::array<Type, size> ConstArray(Type val = 0)
{
    static_assert(size <= std::numeric_limits<std::uint8_t>::max(),
                  "ConstArray loop counter is std::uint8_t and would overflow for this size");
    std::array<Type, size> arr{};
    for (std::uint8_t i = 0; i < size; i++)
        arr[i] = val;
    return arr;
}

template<typename Type>
NO_DISCARD std::vector<Type> displacementFrom(std::vector<Type> const& input)
{
    std::vector<Type> displs(input.size());
    Type off = 0;
    for (Type i = 0; i < static_cast<Type>(input.size()); i++)
    {
        displs[i] = off;
        off += input[i];
    }
    return displs;
}


template<typename F>
NO_DISCARD auto generate_from(F&& f, std::size_t from, std::size_t to)
{
    assert(from <= to);
    using value_type = std::decay_t<std::invoke_result_t<F&, std::size_t const&>>;
    std::vector<value_type> v;
    std::size_t count = to - from;
    if (count > 0)
        v.reserve(count);
    for (std::size_t i = from; i < to; ++i)
        v.emplace_back(f(i));
    return v;
}

template<typename F>
NO_DISCARD auto generate_from(F&& f, std::size_t count)
{
    return generate_from(std::forward<F>(f), 0, count);
}


template<typename F, typename Container>
NO_DISCARD auto generate_from(F&& f, Container const& container)
{
    using T          = typename Container::value_type;
    using value_type = std::decay_t<std::invoke_result_t<F&, T&>>;
    std::vector<value_type> v1;
    if (container.size() > 0)

        v1.reserve(container.size());
    for (auto& v : container)
        v1.emplace_back(f(v));
    return v1;
}

template<typename F, typename T>
NO_DISCARD auto generate_from(F&& f, std::vector<T>&& v)
{
    return generate_from(std::forward<F>(f), v);
}

template<typename F, typename Type, std::size_t Size>
NO_DISCARD auto constexpr generate_from(F&& f, std::array<Type, Size> const& arr)
{
    return for_N_make_array<Size>([&](auto i) { return f(arr[i]); });
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_ARRAY_HPP
