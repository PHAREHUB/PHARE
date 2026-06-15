#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_BASE_PARAMS_HPP
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_BASE_PARAMS_HPP



#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/utilities/point/point.hpp"

#include "core/models/options/hybrid_options.hpp"

#include "gridlayout_utilities.hpp"

#include <memory>

using namespace PHARE::core;

template<auto options>
struct GridLayoutTestParam
{
    PHARE::HybridFieldOptions<options> constexpr static field_options{};
    using GridLayout_t                       = GridLayout<PHARE::HybridOptions<field_options>{}>;
    static constexpr std::size_t dim         = options.dimension;
    static constexpr std::size_t interpOrder = options.interp_order;
    using Grid_t = Grid<decltype(getNdArrayVecImpl(SelectorDim<dim>{})), HybridQuantity::Scalar>;


    std::shared_ptr<GridLayout_t> layout;
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
