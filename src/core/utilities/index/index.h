#ifndef PHARE_CORE_UTILITIES_INDEX_INDEX_H
#define PHARE_CORE_UTILITIES_INDEX_INDEX_H

#include <cstddef>

#include "data/grid/gridlayoutdefs.h"
#include "utilities/types.h"

namespace PHARE
{
template<std::size_t... dim>
struct MeshIndex
{
};


template<>
struct MeshIndex<1>
{
    MeshIndex<1>() = default;

    MeshIndex<1>(uint32 idx, QtyCentering centering)
        : i{idx}
        , centering{centering}
    {
    }
    static constexpr uint32 dimension = 1;
    uint32 i{0};
    QtyCentering centering;
};


template<>
struct MeshIndex<2>
{
    MeshIndex<2>() = default;

    MeshIndex<2>(uint32 idx1, uint32 idx2, QtyCentering centering)
        : i{idx1}
        , j{idx2}
        , centering{centering}
    {
    }
    static constexpr uint32 dimension = 2;
    uint32 i{0}, j{0};
    QtyCentering centering;
};


template<>
struct MeshIndex<3>
{
    MeshIndex<3>() = default;

    MeshIndex<3>(uint32 idx1, uint32 idx2, uint32 idx3, QtyCentering centering)
        : i{idx1}
        , j{idx2}
        , k{idx3}
    {
    }
    static constexpr uint32 dimension = 3;
    uint32 i{0}, j{0}, k{0};
    QtyCentering centering;
};



auto make_index(uint32 i, QtyCentering centering)
{
    return MeshIndex<1>(i, centering);
}

auto make_index(uint32 i, uint32 j, QtyCentering centering)
{
    return MeshIndex<2>(i, j, centering);
}

auto make_index(uint32 i, uint32 j, uint32 k, QtyCentering centering)
{
    return MeshIndex<3>(i, j, k, centering);
}


} // namespace PHARE


#endif
