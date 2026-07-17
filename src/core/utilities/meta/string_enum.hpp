#ifndef PHARE_CORE_UTILITIES_META_STRING_ENUM_HPP
#define PHARE_CORE_UTILITIES_META_STRING_ENUM_HPP

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace PHARE::core
{
// Reflection-lite for scoped enums. If for your new Enum, you locally specialize EnumTraits<Enum>
// with:
//   static constexpr std::string_view label;   // human name, used in error messages
//   static constexpr std::array<std::pair<std::string_view, Enum>, N> names;
//                                               // {lowercase name -> value} parseable options
// then fromString<Enum> / toString<Enum> below work generically.
template<typename Enum>
struct EnumTraits; // primary intentionally left undefined; each enum provides a specialization

template<typename Enum>
Enum fromString(std::string s)
{
    static_assert(std::is_enum_v<Enum>, "fromString requires an enum type");
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    for (auto const& [name, value] : EnumTraits<Enum>::names)
        if (s == name)
            return value;
    throw std::runtime_error("Unknown " + std::string{EnumTraits<Enum>::label} + ": " + s);
}

template<typename Enum>
std::string toString(Enum const value)
{
    static_assert(std::is_enum_v<Enum>, "toString requires an enum type");
    for (auto const& [name, v] : EnumTraits<Enum>::names)
        if (v == value)
            return std::string{name};
    throw std::runtime_error("Unhandled " + std::string{EnumTraits<Enum>::label} + " value");
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_STRING_ENUM_HPP
