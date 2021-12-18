#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_UTILITIES_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_UTILITIES_H

#include <array>

#include "core/data/ndarray/ndarray_vector.hpp"

using namespace PHARE::core;

template<std::size_t dim>
struct SelectorDim
{
};

inline NdArrayVector<1> getNdArrayVecImpl(SelectorDim<1>)
{
    return NdArrayVector<1>{0u};
}

inline NdArrayVector<2> getNdArrayVecImpl(SelectorDim<2>)
{
    return NdArrayVector<2>{0u, 0u};
}

inline NdArrayVector<3> getNdArrayVecImpl(SelectorDim<3>)
{
    return NdArrayVector<3>{0u, 0u, 0u};
}



#endif
