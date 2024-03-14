#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_BASE_PARAMS_HPP
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_BASE_PARAMS_HPP

#include <memory>


#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/point/point.hpp"

#include "gridlayout_utilities.hpp"

using namespace PHARE::core;

template<typename GridLayoutImpl>
struct GridLayoutTestParam
{
    static constexpr std::size_t dim         = GridLayoutImpl::dimension;
    static constexpr std::size_t interpOrder = GridLayoutImpl::interp_order;
    using Grid_t = Grid<decltype(getNdArrayVecImpl(SelectorDim<dim>{})), HybridQuantity::Scalar>;

    std::shared_ptr<GridLayout<GridLayoutImpl>> layout;
    std::array<double, dim> dxdydz;
    std::array<std::uint32_t, dim> nbCellXYZ;

    Point<double, dim> origin;

    HybridQuantity::Scalar currentQuantity;

    std::shared_ptr<Grid_t> field;


    template<typename Container, std::size_t... I>
    auto makeIt_(Container allocSize, std::index_sequence<I...>)
    {
        return std::make_shared<Grid_t>("field", currentQuantity, (allocSize[I])...);
    }

    template<typename Container>
    auto makeMyField_(Container allocSize)
    {
        std::make_index_sequence<dim> idx;
        return makeIt_(allocSize, idx);
    }
};


#endif
