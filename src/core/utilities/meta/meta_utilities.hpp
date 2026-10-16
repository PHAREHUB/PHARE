#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP


#include "core/utilities/types.hpp"

#include <concepts>
#include <stdexcept>
#include <tuple>
#include <cassert>
#include <iterator>
#include <type_traits>
#include <utility>
#include <variant>


namespace PHARE
{
namespace core
{
    template<typename...>
    using tryToInstanciate = void;


    struct dummy
    {
        using type              = int;
        static type const value = 0;
    };


    /** \brief Iterable is satisfied by any type that can be passed to std::begin/std::end,
     * e.g. a Box or a BoxContainer
     */
    template<typename IterableCandidate>
    concept Iterable = requires(IterableCandidate c) {
        std::begin(c);
        std::end(c);
    };


    template<typename IterableCandidate>
    using is_iterable = std::enable_if_t<Iterable<IterableCandidate>, dummy::type>;


    template<typename IterableCandidate>
    constexpr static bool is_iterable_v = Iterable<IterableCandidate>;


    // Basic function
    template<typename T>
    constexpr void allsame(T)
    {
    }

    // Recursive function
    template<typename T, typename T2, typename... Ts,
             typename = std::enable_if_t<std::is_same<T, T2>::value>>
    constexpr void allsame([[maybe_unused]] T arg, T2 arg2, Ts... args)
    {
        allsame(arg2, args...);
    }


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
        return [&]<std::size_t... Is>(std::index_sequence<Is...>) {
            EnumConstantVariant_t<Enum> result;
            bool found = false;
            (void)(((value == static_cast<Enum>(Is))
                        ? (result = std::integral_constant<Enum, static_cast<Enum>(Is)>{},
                           found  = true, true)
                        : false)
                   || ...);
            if (!found)
                throw std::runtime_error("asEnumConstant: enum value out of range");
            return result;
        }(std::make_index_sequence<static_cast<std::size_t>(Enum::count)>{});
    }

    namespace detail
    {
        inline auto toConstexprVariant(std::same_as<bool> auto const value)
        {
            return asBoolConstant(value);
        }

        template<CountedEnum Enum>
        auto toConstexprVariant(Enum const value)
        {
            return asEnumConstant(value);
        }

        template<typename T>
        using ConstexprVariant_t = decltype(toConstexprVariant(std::declval<T const&>()));

        template<typename T>
        constexpr std::size_t nbrConstexprCases = std::variant_size_v<ConstexprVariant_t<T>>;
    } // namespace detail




    /** @brief Lifts a batch of runtime bool/enum values into compile-time constants.
     *
     * Useful when runtime checks are embedded in a compute loop.
     *
     * @note enum arguments must be CountedEnum, i.e. end with a `count` sentinel giving the
     * number of values, so every case can be enumerated.
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
                        detail::toConstexprVariant(vs)...);
                },
                values);
        }

        std::tuple<Args...> values;
    };
    template<typename... Args>
    Constexprifier(Args const&...) -> Constexprifier<Args...>;



} // namespace core

} // namespace PHARE

#endif
