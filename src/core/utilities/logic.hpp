#ifndef PHARE_CORE_UTILITIES_LOGIC_HPP
#define PHARE_CORE_UTILITIES_LOGIC_HPP

#include "core/def.hpp"
#include "core/utilities/meta/type_meta.hpp"

#include <algorithm>

namespace PHARE::core
{
NO_DISCARD auto constexpr any(auto... bools)
    requires(all_are<bool>(bools...))
{
    return (bools || ...);
}

NO_DISCARD auto constexpr all(auto... bools)
    requires(all_are<bool>(bools...))
{
    return (bools && ...);
}

// calls operator bool() or copies bool
auto constexpr static to_bool = [](auto const& v) { return bool{v}; };


template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr all(Container const& container, Fn fn = to_bool)
{
    return std::all_of(container.begin(), container.end(), fn);
}

template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr any(Container const& container, Fn fn = to_bool)
{
    return std::any_of(container.begin(), container.end(), fn);
}


template<typename Container, typename Fn = decltype(to_bool)>
NO_DISCARD auto constexpr none(Container const& container, Fn fn = to_bool)
{
    return std::none_of(container.begin(), container.end(), fn);
}


auto constexpr all_in(auto const a, auto&&... ts)
{
    return ((a == ts) && ...);
}

auto constexpr any_in(auto const a, auto&&... ts)
{
    return ((a == ts) || ...);
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_LOGIC_HPP
