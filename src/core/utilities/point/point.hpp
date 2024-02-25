#ifndef PHARE_CORE_UTILITIES_POINT_POINT_HPP
#define PHARE_CORE_UTILITIES_POINT_POINT_HPP

#include <cassert>
#include <array>
#include <cstddef>
#include <sstream>
#include <ostream>

#include "core/utilities/meta/meta_utilities.hpp"
#include "core/def.hpp"

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
        using value_type                       = Type;

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
            for (std::size_t i = 0; i < dim; ++i)
            {
                r[i] = c[i];
            }
        }

        constexpr Point() { core::fill(Type{0}, r); }

        NO_DISCARD auto& operator[](std::size_t i) { return r[i]; }
        NO_DISCARD auto const& operator[](std::size_t i) const { return r[i]; }


        NO_DISCARD bool operator==(Point const& p) const
        {
            bool areEqual = true;
            for (std::size_t i = 0; i < dim; ++i)
            {
                static_assert(std::is_integral_v<Type>,
                              "this function is only valid for integral type of Point");

                areEqual &= ((*this)[i] == p[i]);
            }
            return areEqual;
        }

        NO_DISCARD bool operator!=(Point const& other) const { return !(*this == other); }


        template<typename DestType = Type>
        NO_DISCARD auto toArray() const
        {
            std::array<DestType, dimension> destArray;
            for (auto i = 0u; i < dimension; ++i)
            {
                destArray[i] = static_cast<DestType>(r[i]);
            }
            return destArray;
        }


        NO_DISCARD std::vector<Type> toVector() const
        {
            return std::vector<Type>(r.data(), r.data() + dimension);
        }

        NO_DISCARD std::string str() const
        {
            std::stringstream ss;
            ss << r[0];
            for (std::size_t i = 1; i < dim; ++i)
            {
                ss << "," << r[i];
            }
            return ss.str();
        }

        NO_DISCARD static Point fromString(std::string csv)
        {
            Point p;
            std::istringstream split(csv);
            std::vector<std::string> tokens;
            for (std::string each; std::getline(split, each, ','); tokens.push_back(each)) {}
            assert(tokens.size() == dimension);
            for (std::size_t i = 0; i < tokens.size(); i++)
            {
                std::stringstream ss;
                ss << tokens[i];
                ss >> p.r[i];
            }
            return p;
        }



        auto operator+(Type const& value) const
        {
            auto copy = *this;
            for (auto iDim = 0u; iDim < dim; ++iDim)
                copy[iDim] += value;
            return copy;
        }
        auto operator+(std::array<Type, dim> const& value) const
        {
            auto copy = *this;
            for (auto iDim = 0u; iDim < dim; ++iDim)
                copy[iDim] += value[iDim];
            return copy;
        }
        auto operator+(Point<Type, dim> const& value) const { return (*this) + value.r; }


        auto operator-(Type const& value) const
        {
            auto copy = *this;
            for (auto iDim = 0u; iDim < dim; ++iDim)
                copy[iDim] -= value;
            return copy;
        }
        auto operator-(std::array<Type, dim> const& value) const
        {
            auto copy = *this;
            for (auto iDim = 0u; iDim < dim; ++iDim)
                copy[iDim] -= value[iDim];
            return copy;
        }
        auto operator-(Point<Type, dim> const& value) const { return (*this) - value.r; }


        NO_DISCARD constexpr auto size() const { return dim; }
        NO_DISCARD auto begin() { return r.begin(); }
        NO_DISCARD auto begin() const { return r.begin(); }
        NO_DISCARD auto end() { return r.end(); }
        NO_DISCARD auto end() const { return r.end(); }

        NO_DISCARD auto& operator*() const { return r; }

    private:
        std::array<Type, dim> r{};
    };

    template<typename... Indexes>
    Point(Indexes... indexes)
        -> Point<typename std::tuple_element<0, std::tuple<Indexes...>>::type, sizeof...(indexes)>;

    template<typename Type, std::size_t dim>
    auto& operator<<(std::ostream& os, Point<Type, dim> const& p)
    {
        os << "( ";
        for (auto& v : p)
            os << v << " ";
        os << ")";
        return os;
    }

} // namespace core
} // namespace PHARE

namespace std
{
template<typename Type, std::size_t dim>
NO_DISCARD PHARE::core::Point<Type, dim> abs(PHARE::core::Point<Type, dim> const& point)
{
    std::array<Type, dim> postive;
    for (std::size_t i = 0; i < dim; ++i)
        postive[i] = std::abs(point[i]);
    return postive;
}

} // namespace std


#endif
