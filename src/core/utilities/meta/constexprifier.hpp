#ifndef PHARE_CORE_UTILITIES_META_CONSTEXPRIFIER_HPP
#define PHARE_CORE_UTILITIES_META_CONSTEXPRIFIER_HPP

#include <tuple>
#include <cstddef>
#include <utility>
#include <variant>
#include <concepts>
#include <type_traits>

namespace PHARE::core
{
/** @brief lifts a runtime bool into a compile-time std::bool_constant (std::false_type /
 * std::true_type) wrapped in a variant, so std::visit can fan out both cases and let the
 * visitor branch via `if constexpr`.
 *
 * argument @c value is declared with concept + auto to avoid implicit conversions to bool.
 */
inline std::variant<std::false_type, std::true_type>
asBoolConstant(std::same_as<bool> auto const value)
{
    if (value)
        return std::true_type{};
    return std::false_type{};
}


/** @brief customization point: how Constexprifier lifts a runtime value of type T into a
 * std::variant of compile-time constants. Specialize it with a static lift(T) to support T;
 * bool is handled here, CountedEnums in enum_meta.hpp. A class template is used, not overloads,
 * so specializations declared after this header are still found when Constexprifier is
 * instantiated.
 */
template<typename T>
struct ConstexprLift; // undefined: T is not liftable

template<>
struct ConstexprLift<bool>
{
    static auto lift(bool const value) { return asBoolConstant(value); }
};

namespace detail
{
    template<typename T>
    using ConstexprVariant_t = decltype(ConstexprLift<T>::lift(std::declval<T const&>()));

    template<typename T>
    constexpr std::size_t nbrConstexprCases = std::variant_size_v<ConstexprVariant_t<T>>;
} // namespace detail




/** @brief Lifts a batch of runtime bool/enum values into compile-time constants.
 *
 * Useful when runtime checks are embedded in a compute loop.
 *
 * @note each argument type needs a ConstexprLift specialization: bool here, CountedEnum (an
 * enum ending with a `count` sentinel, so every case can be enumerated) in enum_meta.hpp.
 *
 * @note the number of visitor instantiations is the *product* of the per-argument case
 * counts, and is capped at MAX_CONSTEXPR_PERMUTATIONS by a static_assert.
 *
 * @code
 * enum class Mode { A, B, C, count };
 *
 * bool aCondition = true;
 * Mode anEnumValue = Mode::B;
 *
 * Constexprifier{aCondition, anEnumValue}([&]<bool constCondition, Mode constEnumValue>() {
 *     for (std::size_t i = 0; i < 10; ++i)
 *     {
 *         if constexpr (constCondition)
 *         {
 *             // do something
 *         }
 *
 *         if constexpr (constEnumValue == Mode::B)
 *         {
 *             // do something
 *         }
 *     }
 * });
 * @endcode
 *
 * @see core::Ohm::operator() in core/numerics/ohm/ohm.hpp (two bools)
 * @see core::Godunov::operator() in core/numerics/godunov_fluxes/godunov_fluxes.hpp
 *      (two bools + an enum)
 */
template<typename... Args>
struct Constexprifier
{
    /** maximum allowed number of compile-time permutations */
    static constexpr std::size_t MAX_CONSTEXPR_PERMUTATIONS = 64;
    /** total number of compile-time instanciations with given arguments */
    static constexpr std::size_t permutations
        = (std::size_t{1} * ... * detail::nbrConstexprCases<Args>);
    static_assert(permutations <= MAX_CONSTEXPR_PERMUTATIONS,
                  "Constexprifier: permutation budget exceeded");

    explicit Constexprifier(Args const&... args)
        : values{args...}
    {
    }

    void operator()(auto&& fn) const
    {
        std::apply(
            [&](auto const&... vs) {
                std::visit(
                    [&](auto... tags) { fn.template operator()<decltype(tags)::value...>(); },
                    ConstexprLift<Args>::lift(vs)...);
            },
            values);
    }

    std::tuple<Args...> values;
};
template<typename... Args>
Constexprifier(Args const&...) -> Constexprifier<Args...>;

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_META_CONSTEXPRIFIER_HPP
