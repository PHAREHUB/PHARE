#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP


#include "core/utilities/types.hpp"
#include "core/utilities/meta/enum.hpp"

#include <concepts>
#include <tuple>
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




    template<typename IterableCandidate, typename AttemptBegin = void, typename AttemptEnd = void>
    struct has_beginend : std::false_type
    {
    };



    /** \brief has_beginend is a traits that permit to check if a Box or a BoxContainer
     * is passed as template argument
     */
    template<typename IterableCandidate>
    struct has_beginend<IterableCandidate,
                        tryToInstanciate<decltype(std::begin(std::declval<IterableCandidate>()))>,
                        tryToInstanciate<decltype(std::end(std::declval<IterableCandidate>()))>>
        : std::true_type
    {
    };


    template<typename IterableCandidate>
    using is_iterable = std::enable_if_t<has_beginend<IterableCandidate>::value, dummy::type>;


    template<typename IterableCandidate>
    constexpr static bool is_iterable_v = has_beginend<IterableCandidate>::value;


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



    namespace detail
    {
        inline auto toConstexprVariant(std::same_as<bool> auto const value)
        {
            return asBoolConstant(value);
        }

        template<typename Enum>
        auto toConstexprVariant(Enum const value)
            requires(std::is_enum_v<Enum>)
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
     * @note enum arguments must have an `EnumTraits<Enum>` specialization (see
     * core/utilities/meta/enum.hpp) listing every value the enum can take.
     *
     * @note the number of visitor instantiations is the *product* of the per-argument case
     * counts, and is capped at MAX_CONSTEXPR_PERMUTATIONS by a static_assert.
     *
     * @code
     * enum class Mode { A, B, C };
     *
     * template<>
     * struct EnumTraits<Mode>
     * {
     *     static constexpr std::string_view label = "Mode";
     *     static constexpr std::array names
     *     {
     *         enumEntry("a", Mode::A),
     *         enumEntry("b", Mode::B),
     *         enumEntry("c", Mode::C)
     *     };
     * };
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


    template<typename DimConstant, typename InterpConstant, std::size_t... ValidNbrParticles>
    using SimulatorOption = std::tuple<DimConstant, InterpConstant,
                                       std::integral_constant<std::size_t, ValidNbrParticles>...>;

    constexpr decltype(auto) possibleSimulators()
    {
        // inner tuple = dim, interp, list[possible nbrParticles for dim/interp]
        return std::tuple<SimulatorOption<DimConst<1>, InterpConst<1>, 2, 3>,
                          SimulatorOption<DimConst<1>, InterpConst<2>, 2, 3, 4>,
                          SimulatorOption<DimConst<1>, InterpConst<3>, 2, 3, 4, 5>,
                          SimulatorOption<DimConst<2>, InterpConst<1>, 4, 5, 8, 9>,
                          SimulatorOption<DimConst<2>, InterpConst<2>, 4, 5, 8, 9, 16>,
                          SimulatorOption<DimConst<2>, InterpConst<3>, 4, 5, 8, 9, 25>,

                          // TODO add in the rest of 3d nbrParticles permutations
                          SimulatorOption<DimConst<3>, InterpConst<1>, 6, 12 /*, 27*/>,
                          SimulatorOption<DimConst<3>, InterpConst<2>, 6, 12>,
                          SimulatorOption<DimConst<3>, InterpConst<3>, 6, 12>

                          >{};
    }


    constexpr std::size_t defaultNbrRefinedParts(std::size_t dim, std::size_t interp)
    {
        auto sims            = possibleSimulators();
        using SimsTuple_t    = decltype(sims);
        auto constexpr nsims = std::tuple_size_v<SimsTuple_t>;

        std::size_t nbRefinedPart = 0;

        for_N<nsims>([&](auto i) {
            using SimOption = std::tuple_element_t<i, SimsTuple_t>;

            if (std::tuple_element_t<0, SimOption>{}() == dim
                and std::tuple_element_t<1, SimOption>{}() == interp)
            {
                nbRefinedPart = std::tuple_element_t<2, SimOption>{};
            }
        });

        assert(nbRefinedPart != 0); // is static_assert in constexpr call

        return nbRefinedPart;
    }



    constexpr decltype(auto) phare_exe_default_simulators()
    {
        // feel free to change as you wish
        return std::tuple<SimulatorOption<DimConst<1>, InterpConst<1>, 2>,
                          SimulatorOption<DimConst<1>, InterpConst<2>, 2>,
                          SimulatorOption<DimConst<1>, InterpConst<3>, 2>,

                          SimulatorOption<DimConst<2>, InterpConst<1>, 4>,
                          SimulatorOption<DimConst<2>, InterpConst<2>, 4>,
                          SimulatorOption<DimConst<2>, InterpConst<3>, 4>,

                          SimulatorOption<DimConst<3>, InterpConst<1>, 6>,
                          SimulatorOption<DimConst<3>, InterpConst<2>, 6>,
                          SimulatorOption<DimConst<3>, InterpConst<3>, 6>

                          >{};
    }


    template<typename Maker> // used from PHARE::amr::Hierarchy
    auto makeAtRuntime(std::size_t dim, Maker&& maker)
    {
        using Ptr_t = decltype(maker(dim, 1));
        Ptr_t p{};

        core::apply(phare_exe_default_simulators(), [&](auto const& simType) {
            using SimuType = std::decay_t<decltype(simType)>;
            using _dim     = typename std::tuple_element<0, SimuType>::type;

            if (!p)
                p = maker(dim, _dim{});
        });

        return p;
    }

    template<typename Maker, typename Pointer, typename Dimension, typename InterpOrder,
             typename... NbRefinedParts>
    void _makeAtRuntime(Maker& maker, Pointer& p, std::size_t userDim, std::size_t userInterpOrder,
                        std::size_t userNbRefinedPart,
                        std::tuple<Dimension, InterpOrder, NbRefinedParts...> const&)
    {
        core::apply(std::tuple<NbRefinedParts...>{}, [&](auto const& nbRefinedPart) {
            if (!p)
                p = maker(userDim, userInterpOrder, userNbRefinedPart, Dimension{}, InterpOrder{},
                          nbRefinedPart);
        });
    }

    template<typename Maker>
    auto makeAtRuntime(std::size_t dim, std::size_t interpOrder, std::size_t nbRefinedPart,
                       Maker&& maker)
    {
        using Ptr_t = decltype(maker(dim, interpOrder, nbRefinedPart, 1, 1, 1));
        Ptr_t p     = nullptr;

        core::apply(phare_exe_default_simulators(), [&](auto const& simType) {
            _makeAtRuntime(maker, p, dim, interpOrder, nbRefinedPart, simType);
        });

        return p;
    }

} // namespace core

} // namespace PHARE

#endif
