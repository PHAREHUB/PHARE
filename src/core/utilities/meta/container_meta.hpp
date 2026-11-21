#ifndef PHARE_CORE_UTILITIES_META_CONTAINER_META_HPP
#define PHARE_CORE_UTILITIES_META_CONTAINER_META_HPP

#include <array>
#include <tuple>
#include <vector>
#include <ranges>
#include <concepts>

namespace PHARE::core
{
template<typename T>
concept IsStdVector
    = std::same_as<T, std::vector<typename T::value_type, typename T::allocator_type>>;

template<typename T> // any size, tuple_size first so other types are false, not errors
concept IsStdArray = requires {
    std::tuple_size<T>::value;
} and std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>>;

/** can be iterated over: begin() and end(), eg a ParticleArray (a Spannable also has data()) */
template<typename T>
concept Iterable = std::ranges::range<T>;

template<typename T, typename Index = int>
concept Subscriptable = requires(T t, Index i) { t[i]; };

/** recursively unwraps ::value_type, eg std::vector<std::array<double, 3>> => double
    Stops: types that have a value_type but should be kept whole rather than unwrapped,
      eg std::string, which would otherwise become char:
        final_type_t<std::vector<std::string>>              => char
        final_type_t<std::vector<std::string>, std::string> => std::string */
template<typename T, typename... Stops>
struct final_type
{
    using type = T;
};

template<typename T, typename... Stops>
    requires(requires { typename T::value_type; } and not(std::same_as<T, Stops> or ...))
struct final_type<T, Stops...>
{
    using type = final_type<typename T::value_type, Stops...>::type;
};

template<typename T, typename... Stops>
using final_type_t = final_type<T, Stops...>::type;

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_CONTAINER_META_HPP
