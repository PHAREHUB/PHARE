#ifndef PHARE_CORE_UTILITIES_POINT_POINT_H
#define PHARE_CORE_UTILITIES_POINT_POINT_H

#include <array>
#include <cstddef>

namespace PHARE
{
template<typename Type, std::size_t dim>
class Point
{
public:
    static constexpr std::size_t dimension = dim;
    using type                             = Type;

    template<typename... Indexes>
    constexpr explicit Point(Indexes... index)
        : r{{index...}}
    {
        static_assert(sizeof...(Indexes) == dimension,
                      "Error dimension does match number of arguments");
    }


    constexpr Point(std::array<Type, dim> coords)
        : r{std::move(coords)}
    {
    }


    constexpr Point() { r.fill(static_cast<Type>(0)); }

    type& operator[](std::size_t i) { return r[i]; }

    type const& operator[](std::size_t i) const { return r[i]; }

private:
    std::array<Type, dim> r;
};




} // namespace PHARE

#endif
