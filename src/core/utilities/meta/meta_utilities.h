#ifndef PHARE_CORE_UTILITIES_META_META_UTILITIES_H
#define PHARE_CORE_UTILITIES_META_META_UTILITIES_H

#include <iterator>
#include <type_traits>

#include "core/utilities/types.h"

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


    template<size_t val>
    using _IC = std::integral_constant<std::size_t, val>;

    static auto constexpr possibleSimulators()
    {
        // inner tuple = dim, interp, list[possible nbrParticles for dim/interp]
        return std::tuple< // formatting
            std::tuple<_IC<1>, _IC<1>, _IC<1>, _IC<2>, _IC<3>>,
            std::tuple<_IC<1>, _IC<2>, _IC<1>, _IC<2>, _IC<3>, _IC<4>>,
            std::tuple<_IC<1>, _IC<3>, _IC<2>, _IC<3>, _IC<4>, _IC<5>>>{};
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
                        size_t userNbRefinedPart,
                        std::tuple<Dimension, InterpOrder, NbRefinedParts...> const&)
    {
        core::apply(std::tuple<NbRefinedParts...>{}, [&](auto const& nbRefinedPart) {
            if (!p)
                p = maker(userDim, userInterpOrder, userNbRefinedPart, Dimension{}, InterpOrder{},
                          nbRefinedPart);
        });
    }

    template<typename Maker>
    auto makeAtRuntime(std::size_t dim, std::size_t interpOrder, size_t nbRefinedPart,
                       Maker&& maker)
    {
        using Ptr_t = decltype(maker(dim, interpOrder, nbRefinedPart, 1, 1, 1));
        Ptr_t p{};

        core::apply(possibleSimulators(), [&](auto const& simType) {
            _makeAtRuntime(maker, p, dim, interpOrder, nbRefinedPart, simType);
        });

        return p;
    }

} // namespace core

} // namespace PHARE

#endif
