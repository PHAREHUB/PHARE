#ifndef PHARE_CORE_UTILITIES_META_TYPE_META_HPP
#define PHARE_CORE_UTILITIES_META_TYPE_META_HPP

#include <tuple>
#include <cstddef>
#include <utility>
#include <type_traits>

namespace PHARE::core
{
template<typename>
inline constexpr bool dependent_false_v = false;


template<typename T, std::size_t... Is>
constexpr auto gft_helper(std::index_sequence<Is...> const&&)
    -> decltype(std::make_tuple((Is, std::declval<T>())...));

template<typename T, std::size_t N>
constexpr auto get_fixed_tuple() -> decltype(gft_helper<T>(std::make_index_sequence<N>{}));

template<typename T, std::size_t N>
using tuple_fixed_type = decltype(get_fixed_tuple<T, N>());


template<typename T>
auto constexpr all_are(auto&&... ts)
{
    return ((std::is_same_v<T, std::decay_t<decltype(ts)>>) && ...);
}


template<typename... Ts>
auto constexpr is_any_of(auto const& t)
{
    return ((std::is_same_v<Ts, std::decay_t<decltype(t)>>) || ...);
}

template<typename T, typename... Ts>
auto constexpr is_any_of()
{
    return ((std::is_same_v<Ts, T>) || ...);
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_TYPE_META_HPP
