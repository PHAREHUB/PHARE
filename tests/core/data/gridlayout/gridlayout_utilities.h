#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_UTILITIES_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_UTILITIES_H

#include <array>

#include "core/data/ndarray/ndarray_vector.h"

using namespace PHARE::core;

template<std::size_t dim>
struct SelectorDim
{
};

inline NdArrayVector1D<> getNdArrayVecImpl(SelectorDim<1>)
{
    return NdArrayVector1D<>{0};
}

inline NdArrayVector2D<> getNdArrayVecImpl(SelectorDim<2>)
{
    return NdArrayVector2D<>{0, 0};
}

inline NdArrayVector3D<> getNdArrayVecImpl(SelectorDim<3>)
{
    return NdArrayVector3D<>{0, 0, 0};
}



#endif
