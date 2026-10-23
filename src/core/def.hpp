#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <cassert>
#include <utility>
#include <stdexcept>
#include <type_traits>

#include "core/def/phare_config.hpp" // IWYU pragma: keep

#define NO_DISCARD [[nodiscard]]

#if !defined(NDEBUG) || defined(PHARE_FORCE_DEBUG_DO)
#define PHARE_DEBUG_DO(...) __VA_ARGS__
#define PHARE_DEBUG 1
#else
#define PHARE_DEBUG_DO(...)
#define PHARE_DEBUG 0
#endif


#if !defined(PHARE_UNDEF_ASSERT)
//  Cuda can fail to compile with assertions
//  I've seen a github issue, will ref
#define PHARE_ASSERT(...) assert(__VA_ARGS__)
#else
#define PHARE_ASSERT(...)
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
    static_assert(sizeof...(args) > 0);
    auto constexpr check = [](auto const& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            return arg != nullptr;
        else
            return arg.isUsable();
    };
    return (check(args) && ...);
}


NO_DISCARD bool isSettable(auto const&... args)
{
    static_assert(sizeof...(args) > 0);
    auto constexpr check = [](auto const& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<decltype(arg)>>)
            return arg == nullptr;
        else
            return arg.isSettable();
    };
    return (check(args) && ...);
}

} // namespace PHARE::core


namespace PHARE
{
// device code can't throw: args are unused there and it asserts instead
template<typename Exception, typename... Args>
inline void throw_exception([[maybe_unused]] Args&&... args) _PHARE_ALL_FN_
{
#if defined(__HIPCC__) || defined(__CUDACC__)
    PHARE_ASSERT(false);
#else
    throw Exception{std::forward<Args>(args)...};
#endif
}

template<typename T>
inline void throw_runtime_error(T const& err) _PHARE_ALL_FN_
{
    throw_exception<std::runtime_error>(err);
}

} // namespace PHARE



#endif // PHARE_CORE_DEF_HPP
