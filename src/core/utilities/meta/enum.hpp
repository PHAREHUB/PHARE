#ifndef PHARE_CORE_UTILITIES_META_ENUM_HPP
#define PHARE_CORE_UTILITIES_META_ENUM_HPP

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace PHARE::core
{
/** @brief reflection-lite for scoped enums: locally specialize EnumTraits<Enum> with
 *   static constexpr std::string_view label;   // human name, used in error messages
 *   static constexpr std::array<std::pair<std::string_view, Enum>, N> names;
 *                                               // {lowercase name -> value} parseable options
 * then fromString<Enum> / toString<Enum> / asEnumConstant<Enum> below work generically.
 *
 * The primary template is intentionally left undefined; each enum provides a specialization.
 */
template<typename Enum>
struct EnumTraits;

/** @brief parses a (case-insensitive) string into its EnumTraits<Enum>::names value, throwing
 * if no name matches.
 */
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

/** @brief returns the EnumTraits<Enum>::names name of a value, throwing if unlisted. */
template<typename Enum>
std::string toString(Enum const value)
{
    static_assert(std::is_enum_v<Enum>, "toString requires an enum type");
    for (auto const& [name, v] : EnumTraits<Enum>::names)
        if (v == value)
            return std::string{name};
    throw std::runtime_error("Unhandled " + std::string{EnumTraits<Enum>::label} + " value");
}


/** @brief undefined helper whose return type is the variant of std::integral_constant over the
 * values listed in EnumTraits<Enum>::names; used only to spell EnumConstantVariant_t.
 */
template<typename Enum, std::size_t... Is>
auto enumConstantVariant(std::index_sequence<Is...>)
    -> std::variant<std::integral_constant<Enum, EnumTraits<Enum>::names[Is].second>...>;

/** @brief the std::variant of compile-time tags returned by asEnumConstant<Enum>. */
template<typename Enum>
using EnumConstantVariant_t = decltype(enumConstantVariant<Enum>(
    std::make_index_sequence<EnumTraits<Enum>::names.size()>{}));


/** @brief implementation of asEnumConstant: scans EnumTraits<Enum>::names for the matching
 * value and returns the corresponding variant alternative, throwing if none matches.
 */
template<typename Enum, std::size_t... Is>
EnumConstantVariant_t<Enum> asEnumConstant_(Enum const value, std::index_sequence<Is...>)
{
    EnumConstantVariant_t<Enum> result;
    bool found = false;
    (void)(((value == EnumTraits<Enum>::names[Is].second)
                ? (result = std::integral_constant<Enum, EnumTraits<Enum>::names[Is].second>{},
                   found  = true, true)
                : false)
           || ...);
    if (!found)
        throw std::runtime_error("Unknown " + std::string{EnumTraits<Enum>::label});
    return result;
}


/** @brief lifts a runtime enum value into a compile-time std::integral_constant wrapped in a
 * variant over EnumTraits<Enum>::names, so std::visit can fan out every case and let the
 * visitor branch via `if constexpr`. Values and label come from EnumTraits<Enum> — no set to
 * restate at the call site.
 *
 * Requires a defined EnumTraits<Enum> specialization; without one this fails to compile.
 */
template<typename Enum>
auto asEnumConstant(Enum const value)
{
    static_assert(std::is_enum_v<Enum>, "asEnumConstant requires an enum type");
    return asEnumConstant_<Enum>(value,
                                 std::make_index_sequence<EnumTraits<Enum>::names.size()>{});
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_ENUM_HPP
