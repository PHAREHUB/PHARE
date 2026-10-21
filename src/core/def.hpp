#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <type_traits>

#define NO_DISCARD [[nodiscard]]

// for helpers called inside loops that must vectorise, where the inliner's cost model otherwise
// keeps them out of line: PHARE_LAMBDA_ALWAYS_INLINE goes after a lambda's parameters,
// PHARE_ALWAYS_INLINE before a function's return type
#if defined(__GNUC__) || defined(__clang__)
#define PHARE_LAMBDA_ALWAYS_INLINE __attribute__((always_inline))
#define PHARE_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define PHARE_LAMBDA_ALWAYS_INLINE
#define PHARE_ALWAYS_INLINE inline
#endif

// on a loop whose iterations are independent but which reads and writes more arrays than the
// compiler checks for overlap at run time: lets it vectorise without those checks
#if defined(__clang__)
#define PHARE_IVDEP _Pragma("clang loop vectorize(assume_safety)")
#elif defined(__GNUC__)
#define PHARE_IVDEP _Pragma("GCC ivdep")
#else
#define PHARE_IVDEP
#endif

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
