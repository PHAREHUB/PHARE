#ifndef PHARE_CORE_UTILITIES_INDEX_INDEX_H
#define PHARE_CORE_UTILITIES_INDEX_INDEX_H

#include <cstddef>

#include "data/grid/gridlayoutdefs.h"
#include "utilities/types.h"

namespace PHARE
{
namespace core


{
    template<std::size_t... dim>
    struct MeshIndex
    {
    };


    template<>
    struct MeshIndex<1>
    {
        MeshIndex<1>() = default;

        MeshIndex<1>(uint32 idx)
            : i{idx}
        {
        }
        static constexpr uint32 dimension = 1;
        uint32 i{0};
    };


    template<>
    struct MeshIndex<2>
    {
        MeshIndex<2>() = default;

        MeshIndex<2>(uint32 idx1, uint32 idx2)
            : i{idx1}
            , j{idx2}
        {
        }
        static constexpr uint32 dimension = 2;
        uint32 i{0}, j{0};
    };


    template<>
    struct MeshIndex<3>
    {
        MeshIndex<3>() = default;

        MeshIndex<3>(uint32 idx1, uint32 idx2, uint32 idx3)
            : i{idx1}
            , j{idx2}
            , k{idx3}
        {
        }
        static constexpr uint32 dimension = 3;
        uint32 i{0}, j{0}, k{0};
    };



    MeshIndex<1> make_index(uint32 i);
    MeshIndex<2> make_index(uint32 i, uint32 j);
    MeshIndex<3> make_index(uint32 i, uint32 j, uint32 k);


} // namespace core
} // namespace PHARE


#endif
