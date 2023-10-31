#ifndef PHARE_CORE_UTILITIES_INDEX_INDEX_HPP
#define PHARE_CORE_UTILITIES_INDEX_INDEX_HPP

#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"

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
    [[nodiscard]] MeshIndex<1> make_index(std::uint32_t i);
    [[nodiscard]] MeshIndex<2> make_index(std::uint32_t i, std::uint32_t j);
    [[nodiscard]] MeshIndex<3> make_index(std::uint32_t i, std::uint32_t j, std::uint32_t k);


} // namespace core
} // namespace PHARE


#endif
