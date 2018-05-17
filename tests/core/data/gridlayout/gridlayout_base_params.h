#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_BASE_PARAMS_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_BASE_PARAMS_H

#include <memory>


#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutdefs.h"
#include "data/ndarray/ndarray_vector.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"

namespace PHARE
{
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
struct GridLayoutTestParam
{
    std::shared_ptr<GridLayout<GridLayoutImpl, dim, interpOrder>> layout;

    std::array<double, dim> dxdydz;
    std::array<uint32, dim> nbCellXYZ;

    Point<double, dim> origin;

    HybridQuantity::Scalar currentQuantity;

    std::shared_ptr<Field<decltype(getNdArrayVecImpl(SelectorDim<dim>{})), HybridQuantity::Scalar>>
        field;


    template<typename Container, std::size_t... I>
    auto makeIt_(Container allocSize, std::index_sequence<I...>)
    {
        return std::make_shared<
            Field<decltype(getNdArrayVecImpl(SelectorDim<dim>{})), HybridQuantity::Scalar>>(
            "field", currentQuantity, (allocSize[I])...);
    }

    template<typename Container>
    auto makeMyField_(Container allocSize)
    {
        std::make_index_sequence<dim> idx;
        return makeIt_(allocSize, idx);
    }
};

} // namespace PHARE
#endif
