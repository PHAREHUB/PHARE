#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_HPP

#include <iterator>
#include <type_traits>

#include "core/utilities/types.hpp"

#if !defined(PHARE_SIMULATORS)
#define PHARE_SIMULATORS 3
#endif

namespace PHARE
{
namespace core
{
    template<typename...>
    using tryToInstanciate = void;


    struct dummy
    {
        using type              = int;
        static const type value = 0;
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


    template<typename DimConstant, typename InterpConstant, std::size_t... ValidNbrParticles>
    using SimulatorOption = std::tuple<DimConstant, InterpConstant,
                                       std::integral_constant<std::size_t, ValidNbrParticles>...>;

    constexpr decltype(auto) possibleSimulators()
    {
        // inner tuple = dim, interp, list[possible nbrParticles for dim/interp]
        return std::tuple<SimulatorOption<DimConst<1>, InterpConst<1>, 2, 3>,
                          SimulatorOption<DimConst<1>, InterpConst<2>, 2, 3, 4>,
                          SimulatorOption<DimConst<1>, InterpConst<3>, 2, 3, 4, 5>

#if PHARE_SIMULATORS > 1
                          ,
                          SimulatorOption<DimConst<2>, InterpConst<1>, 4, 5, 8, 9>,
                          SimulatorOption<DimConst<2>, InterpConst<2>, 4, 5, 8, 9, 16>,
                          SimulatorOption<DimConst<2>, InterpConst<3>, 4, 5, 8, 9, 25>
#endif

                          // TODO add in the rest of 3d nbrParticles permutations
                          // possibly consider compile time activation for uncommon cases
#if PHARE_SIMULATORS > 2
                          ,
                          SimulatorOption<DimConst<3>, InterpConst<1>, 6, 12 /*, 27*/>,
                          SimulatorOption<DimConst<3>, InterpConst<2>, 6, 12>,
                          SimulatorOption<DimConst<3>, InterpConst<3>, 6, 12>
#endif
>{};
    }


    template<typename Maker> // used from PHARE::amr::Hierarchy
    auto makeAtRuntime(std::size_t dim, Maker&& maker)
    {
        using Ptr_t = decltype(maker(dim, 1));
        Ptr_t p{};

        core::apply(possibleSimulators(), [&](auto const& simType) {
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

        core::apply(possibleSimulators(), [&](auto const& simType) {
            _makeAtRuntime(maker, p, dim, interpOrder, nbRefinedPart, simType);
        });

        return p;
    }

} // namespace core

} // namespace PHARE

#endif
