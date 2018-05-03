#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_UTILITIES_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_UTILITIES_H

#include <array>

#include "data/ndarray/ndarray_vector.h"

namespace PHARE
{
template<std::size_t dim>
struct SelectorDim
{
};

inline NdArrayVector1D<> getNdArrayVecImpl(SelectorDim<1> dim)
{
    return NdArrayVector1D<>{0};
}

inline NdArrayVector2D<> getNdArrayVecImpl(SelectorDim<2> dim)
{
    return NdArrayVector2D<>{0, 0};
}

inline NdArrayVector3D<> getNdArrayVecImpl(SelectorDim<3> dim)
{
    return NdArrayVector3D<>{0, 0, 0};
}

} // namespace PHARE

#endif
