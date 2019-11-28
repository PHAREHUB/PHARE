#ifndef PHARE_CORE_UTILITIES_POINT_POINT_H
#define PHARE_CORE_UTILITIES_POINT_POINT_H

#include <array>
#include <cstddef>

#include "utilities/meta/meta_utilities.h"

namespace PHARE
{
namespace core
{
    template<typename T, typename Index, typename Attempt = void>
    struct has_subscript_operator : std::false_type
    {
    };


    template<typename T, typename Index>
    struct has_subscript_operator<
        T, Index, tryToInstanciate<decltype(std::declval<T>()[std::declval<Index>()])>>
        : std::true_type
    {
    };


    template<typename T, typename Index = int>
    using is_subscriptable = std::enable_if_t<has_subscript_operator<T, Index>::value, dummy::type>;


    template<typename Type, std::size_t dim>
    class Point
    {
    public:
        static constexpr std::size_t dimension = dim;
        using type                             = Type;

        template<typename... Indexes>
        constexpr Point(Indexes... index)
            : r{{index...}}
        {
            allsame(index...);
            static_assert(sizeof...(Indexes) == dimension,
                          "Error dimension does match number of arguments");
        }


        constexpr Point(std::array<Type, dim> coords)
            : r{std::move(coords)}
        {
        }

        template<typename Container, is_subscriptable<Container> = dummy::value>
        Point(Container c)
        {
            for (auto i = 0u; i < dim; ++i)
            {
                r[i] = c[i];
            }
        }

        constexpr Point() { r.fill(static_cast<Type>(0)); }

        type& operator[](std::size_t i) { return r[i]; }

        type const& operator[](std::size_t i) const { return r[i]; }


        bool operator==(Point const& p) const
        {
            bool areEqual = true;
            for (auto i = 0u; i < dim; ++i)
            {
                static_assert(std::is_integral_v<Type>,
                              "this function is only valid for integral type of Point");

                areEqual &= ((*this)[i] == p[i]);
            }
            return areEqual;
        }


        template<typename DestType>
        auto toArray() const
        {
            std::array<DestType, dimension> destArray;
            for (auto i = 0u; i < dimension; ++i)
            {
                destArray[i] = static_cast<DestType>(r[i]);
            }
            return destArray;
        }



    private:
        std::array<Type, dim> r;
    };

    template<typename... Indexes>
    Point(Indexes... indexes)
        ->Point<typename std::tuple_element<0, std::tuple<Indexes...>>::type, sizeof...(indexes)>;


} // namespace core
} // namespace PHARE

#endif
