#ifndef PHARE_CORE_UTILITIES_META_ENUM_META_HPP
#define PHARE_CORE_UTILITIES_META_ENUM_META_HPP

#include "core/utilities/meta/for_N.hpp"
#include "core/utilities/meta/constexprifier.hpp"

#include <cstdint>
#include <utility>
#include <variant>
#include <stdexcept>
#include <type_traits>

namespace PHARE::core
{
template<typename Enum>
concept CountedEnum = std::is_enum_v<Enum> and requires { Enum::count; };

/**
 * @brief undefined helper whose return type is the variant of std::integral_constant over
 * every value of a CountedEnum; used only to spell EnumConstantVariant_t.
 */
template<typename Enum, std::size_t... Is>
auto enumConstantVariant(std::index_sequence<Is...>)
    -> std::variant<std::integral_constant<Enum, static_cast<Enum>(Is)>...>;

/** @brief alias for the std::variant of compile-time tags returned by asEnumConstant<Enum>. */
template<CountedEnum Enum>
using EnumConstantVariant_t = decltype(enumConstantVariant<Enum>(
    std::make_index_sequence<static_cast<std::size_t>(Enum::count)>{}));

/** @brief transforms a runtime enum value into a compile-time std::integral_constant wrapped in
 * a variant over every value of the enum. Throws if the value is not below Enum::count.
 */
template<CountedEnum Enum>
EnumConstantVariant_t<Enum> asEnumConstant(Enum const value)
{
    return for_N_fold<static_cast<std::uint16_t>(Enum::count)>([&](auto... is) {
        EnumConstantVariant_t<Enum> result;
        bool found = false;
        (void)(((value == static_cast<Enum>(is()))
                    ? (result = std::integral_constant<Enum, static_cast<Enum>(is())>{},
                       found  = true, true)
                    : false)
               || ...);
        if (!found)
            throw std::runtime_error("asEnumConstant: enum value out of range");
        return result;
    });
}

/** lets Constexprifier lift a CountedEnum, see asEnumConstant */
template<CountedEnum Enum>
struct ConstexprLift<Enum>
{
    static auto lift(Enum const value) { return asEnumConstant(value); }
};

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_ENUM_META_HPP
