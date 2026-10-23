#ifndef PHARE_CORE_UTILITIES_INDEX_INDEX_HPP
#define PHARE_CORE_UTILITIES_INDEX_INDEX_HPP

#include <array>
#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"
#include "core/def.hpp"

namespace PHARE
{
namespace core
{
    template<std::size_t dim>
    using MeshIndex = Point<std::size_t, dim>;



    /*
        template<std::size_t... dim>
        struct MeshIndex
        {
        };

            template<>
            struct MeshIndex<1>
            {
                MeshIndex<1>() = default;

                MeshIndex<1>(std::uint32_t idx)
                    : i{idx}
                {
                }
                static constexpr std::uint32_t dimension = 1;
                std::uint32_t i{0};
            };


            template<>
            struct MeshIndex<2>
            {
                MeshIndex<2>() = default;

                MeshIndex<2>(std::uint32_t idx1, std::uint32_t idx2)
                    : i{idx1}
                    , j{idx2}
                {
                }
                static constexpr std::uint32_t dimension = 2;
                std::uint32_t i{0}, j{0};
            };


            template<>
            struct MeshIndex<3>
            {
                MeshIndex<3>() = default;

                MeshIndex<3>(std::uint32_t idx1, std::uint32_t idx2, std::uint32_t idx3)
                    : i{idx1}
                    , j{idx2}
                    , k{idx3}
                {
                }
                static constexpr std::uint32_t dimension = 3;
                std::uint32_t i{0}, j{0}, k{0};
            };


        */
    // inline: stencil neighbours (GridLayout::next/previous) build one of these per access
    NO_DISCARD inline MeshIndex<1> make_index(std::uint32_t i)
    {
        return MeshIndex<1>(i);
    }

    NO_DISCARD inline MeshIndex<2> make_index(std::uint32_t i, std::uint32_t j)
    {
        return MeshIndex<2>(i, j);
    }

    NO_DISCARD inline MeshIndex<3> make_index(std::uint32_t i, std::uint32_t j, std::uint32_t k)
    {
        return MeshIndex<3>(i, j, k);
    }



    // calls fn(MeshIndex) for every index of the inclusive box [lo, hi], last dimension fastest
    template<std::size_t dim>
    void for_each_index(std::array<int, dim> const& lo, std::array<int, dim> const& hi, auto&& fn)
    {
        for (int i = lo[0]; i <= hi[0]; ++i)
            if constexpr (dim == 1)
                fn(make_index(i));
            else
                for (int j = lo[1]; j <= hi[1]; ++j)
                    if constexpr (dim == 2)
                        fn(make_index(i, j));
                    else
                        for (int k = lo[2]; k <= hi[2]; ++k)
                            fn(make_index(i, j, k));
    }

} // namespace core
} // namespace PHARE


#endif
