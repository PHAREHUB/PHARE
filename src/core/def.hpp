#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <type_traits>

#define NO_DISCARD [[nodiscard]]

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__
#else
#define PHARE_DEBUG_DO(...)
#endif

#define _PHARE_TO_STR(x) #x // convert macro text to string
#define PHARE_TO_STR(x) _PHARE_TO_STR(x)

#define PHARE_TOKEN_PASTE(x, y) x##y
#define PHARE_STR_CAT(x, y) PHARE_TOKEN_PASTE(x, y)


namespace PHARE::core
{


template<typename T>
concept FloatingPoint = std::is_floating_point_v<T>;



NO_DISCARD bool isUsable(auto const&... args)
{
    auto check = [](auto const& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            return arg != nullptr;
        else
            return arg.isUsable();
    };
    return (check(args) && ...);
}


NO_DISCARD bool isSettable(auto const&... args)
{
    auto check = [](auto const& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            return arg == nullptr;
        else
            return arg.isSettable();
    };
    return (check(args) && ...);
}

} // namespace PHARE::core

#endif // PHARE_CORE_DEF_HPP
