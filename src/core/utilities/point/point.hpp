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


        constexpr Point(std::array<Type, dim> const& coords)
            : r{coords}
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


        template<typename T2>
        NO_DISCARD bool operator==(Point<T2, dim> const& p) const
        {
            bool areEqual = true;
            for (std::size_t i = 0; i < dim; ++i)
            {
                static_assert(std::is_integral_v<Type>,
                              "this function is only valid for integral type of Point");

                areEqual &= int_equals((*this)[i], p[i]); // handles signed differences
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

        NO_DISCARD std::string str(std::uint8_t const precision = 6,
                                   bool left_pad_positive       = true) const
            requires(FloatingPoint<Type>)
        {
            std::stringstream ss;
            for (auto const v : r)
            {
                if (left_pad_positive and v >= 0)
                    ss << " ";
                ss << to_string_with_precision(v, precision) << ",";
            }
            auto s = ss.str(); // drop last comma
            s.pop_back();
            return s;
        }


        NO_DISCARD std::string str() const
            requires(!FloatingPoint<Type>)
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

        auto& operator+=(Type const& value)
        {
            for (auto iDim = 0u; iDim < dim; ++iDim)
                r[iDim] += value;
            return *this;
        }

        template<template<typename, std::size_t> typename Arr, typename T>
        auto& operator+=(Arr<T, dim> const& value)
        {
            for (auto iDim = 0u; iDim < dim; ++iDim)
                r[iDim] += value[iDim];
            return *this;
        }
        auto& operator-=(Type const& value)
        {
            for (auto iDim = 0u; iDim < dim; ++iDim)
                r[iDim] -= value;
            return *this;
        }
        template<template<typename, std::size_t> typename Arr, typename T>
        auto& operator-=(Arr<T, dim> const& value)
        {
            for (auto iDim = 0u; iDim < dim; ++iDim)
                r[iDim] -= value[iDim];
            return *this;
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

        auto operator*(Type const& value) const
        {
            auto copy = *this;
            for (auto iDim = 0u; iDim < dim; ++iDim)
                copy[iDim] *= value;
            return copy;
        }
        auto operator*(std::array<Type, dim> const& value) const
        {
            auto copy = *this;
            for (auto iDim = 0u; iDim < dim; ++iDim)
                copy[iDim] *= value[iDim];
            return copy;
        }
        auto operator*(Point<Type, dim> const& value) const { return (*this) * value.r; }


        NO_DISCARD constexpr auto size() const { return dim; }
        NO_DISCARD auto begin() { return r.begin(); }
        NO_DISCARD auto begin() const { return r.begin(); }
        NO_DISCARD auto end() { return r.end(); }
        NO_DISCARD auto end() const { return r.end(); }

        NO_DISCARD auto& operator*() const { return r; }



        template<typename To>
        auto as() const
        {
            return Point<To, dim>{this->template toArray<To>()};
        }

        auto as_unsigned() const
        {
            for (auto iDim = 0u; iDim < dim; ++iDim)
                if (r[iDim] < 0)
                    throw std::runtime_error("Cannot make unsigned from negative values");
            if constexpr (sizeof(Type) == 4)
                return as<std::uint32_t>();
            // else no return cause not yet handled
        }

        auto as_signed() const
        {
            if constexpr (sizeof(Type) == 4)
                return as<std::int32_t>();
            // else no return cause not yet handled
        }

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
