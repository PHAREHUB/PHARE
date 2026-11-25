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

/** the innermost value type of nested std::vector / std::array, eg
      std::vector<std::array<double, 3>> => double
    other types are kept whole: std::vector<std::string> => std::string, std::optional<int> */
template<typename T>
struct final_type
{
    using type = T;
};

template<typename T>
    requires(IsStdVector<T> or IsStdArray<T>)
struct final_type<T>
{
    using type = final_type<typename T::value_type>::type;
};

template<typename T>
using final_type_t = final_type<T>::type;

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_CONTAINER_META_HPP
