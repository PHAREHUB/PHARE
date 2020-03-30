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




    using OneD   = std::integral_constant<std::size_t, 1>;
    using TwoD   = std::integral_constant<std::size_t, 2>;
    using ThreeD = std::integral_constant<std::size_t, 3>;

    using FirstOrder  = std::integral_constant<std::size_t, 1>;
    using SecondOrder = std::integral_constant<std::size_t, 2>;
    using ThirdOrder  = std::integral_constant<std::size_t, 3>;




    static auto constexpr possibleDimensions() { return std::tuple<OneD, TwoD, ThreeD>{}; }


    static auto constexpr possibleInterpOrders()
    {
        return std::tuple<FirstOrder, SecondOrder, ThirdOrder>{};
    }

    template<size_t dim, size_t interpOrder>
    static auto constexpr possibleNmbSplitPartiles()
    {
        if constexpr (dim == 1)
        {
            return std::tuple<OneD, TwoD, ThreeD>{};
        }
        return std::tuple<OneD, TwoD, ThreeD>{};
    }

    template<typename Maker>
    auto makeAtRuntime(std::size_t dim, std::size_t interpOrder, Maker&& maker)
    {
        using Ptr_t = decltype(maker(dim, interpOrder, 1, 1));
        Ptr_t p;
        if constexpr (std::is_same_v<Ptr_t, bool>)
            p = false;
        else
            p = nullptr;

        core::apply(possibleDimensions(), [&](auto const& dim1) {
            core::apply(possibleInterpOrders(), [&](auto const& order) {
                if (!p)
                    p = maker(dim, interpOrder, dim1, order);
            });
        });

        return p;
    }

    template<typename Maker, typename Pointer, typename Dimension, typename InterpOrder>
    auto _makeAtRuntime(Maker& maker, Pointer& p, std::size_t userDim, std::size_t userInterpOrder,
                        size_t userNbRefinedPart, Dimension dimension, InterpOrder interp_order)
    {
        constexpr size_t dim    = dimension();
        constexpr size_t interp = interp_order();

        constexpr auto nmbSplitPartiles = possibleNmbSplitPartiles<dim, interp>();
        core::apply(nmbSplitPartiles, [&](auto const& nbRefinedPart) {
            if (!p)
                p = maker(userDim, userInterpOrder, userNbRefinedPart, dimension, interp_order,
                          nbRefinedPart);
        });
    }


    template<typename Maker>
    auto makeAtRuntime(std::size_t dim, std::size_t interpOrder, size_t nbRefinedPart,
                       Maker&& maker)
    {
        using Ptr_t = decltype(maker(dim, interpOrder, nbRefinedPart, 1, 1, 1));
        Ptr_t p;
        if constexpr (std::is_same_v<Ptr_t, bool>)
            p = false;
        else
            p = nullptr;

        core::apply(possibleDimensions(), [&](auto const& dim1) {
            core::apply(possibleInterpOrders(), [&](auto const& order) {
                _makeAtRuntime(maker, p, dim, interpOrder, nbRefinedPart, dim1, order);
            });
        });

        return p;
    }


} // namespace core

} // namespace PHARE




#endif
