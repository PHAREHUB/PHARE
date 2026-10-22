#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <string>
#include <iostream>
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
concept HasNameMethod = requires(Resource const& res) { res.name(); };

template<typename Resource>
concept HasArrowNameMethod = requires(Resource const& res) { res->name(); };

auto get_resource_name(auto const& res)
    requires HasNameMethod<std::decay_t<decltype(res)>>
{
    return res.name();
}

auto get_resource_name(auto const& res)
    requires(!HasNameMethod<std::decay_t<decltype(res)>>
             and HasArrowNameMethod<std::decay_t<decltype(res)>>)
{
    return res->name();
}

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
                std::cerr << __FILE__ << ":" << __LINE__ << " - "
                          << detail::get_resource_name(arg) << " not usable!" << std::endl;
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
                std::cerr << __FILE__ << ":" << __LINE__ << " - "
                          << detail::get_resource_name(arg) << " not settable!" << std::endl;
        })
        return settable;
    };
    return (check(args) && ...);
}

} // namespace PHARE::core

#endif // PHARE_CORE_DEF_HPP
