#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include "core/logger.hpp"


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



namespace PHARE::core::detail
{
template<typename Resource>
auto get_resource_name(Resource const& res) -> decltype(std::declval<Resource>().name())
{
    return res.name();
}

template<typename Resource>
auto get_resource_name(Resource const* res) -> decltype(std::declval<Resource>()->name())
{
    return res->name();
}

template<typename... Args>
auto get_resource_name(auto const&...)
{
    return std::string{"unknown resource"};
}

} // namespace PHARE::core::detail


namespace PHARE::core
{


template<typename T>
concept FloatingPoint = std::is_floating_point_v<T>;



NO_DISCARD bool isUsable(auto const&... args)
{
    auto check = [](auto const& arg) {
        bool usable = true;
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            usable = arg != nullptr;
        else
            usable = arg.isUsable();
        PHARE_DEBUG_DO({
            if (!usable)
                PHARE_LOG_LINE_SS(detail::get_resource_name(arg) << " not usable!");
        })
        return usable;
    };
    return (check(args) && ...);
}


NO_DISCARD bool isSettable(auto const&... args)
{
    auto check = [](auto const& arg) {
        bool settable = true;
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            settable = arg == nullptr;
        else
            settable = arg.isSettable();
        PHARE_DEBUG_DO({
            if (!settable)
                PHARE_LOG_LINE_SS(detail::get_resource_name(arg) << " not settable!");
        })
        return settable;
    };
    return (check(args) && ...);
}

} // namespace PHARE::core

#endif // PHARE_CORE_DEF_HPP
