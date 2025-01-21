#ifndef PHARE_CORE_DEF_HPP
#define PHARE_CORE_DEF_HPP

#include <bitset>
#include <cstdint>
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


#ifndef PHARE_DOUBLE_BITSET
constexpr std::bitset<5> doubles{0b11110}; // index 0 starts on the right in binary
#else                                      // PHARE_DOUBLE_BITSET
constexpr std::bitset<5> doubles{PHARE_DOUBLE_BITSET};
#endif

template<std::uint8_t i>
bool constexpr is_double()
{
    // 0 = particle delta
    // 1 = particle v
    // 2 = particle charge
    // 3 = particle weight
    // 4 = fields

    return doubles[i] == 1;
}

template<std::uint8_t i>
struct Floater
{
    using value_type = std::conditional_t<is_double<i>(), double, float>;
};

template<std::uint8_t i>
using floater_t = Floater<i>::value_type;


} // namespace PHARE::core

#endif // PHARE_CORE_DEF_HPP
