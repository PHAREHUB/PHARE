
#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_HPP

#include <array>

#include "core/data/grid/gridlayout.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/utilities/point/point.hpp"

using namespace PHARE::core;


template<typename GridLayoutImpl>
struct GridLayoutCellCenteringParam
{
    GridLayoutTestParam<GridLayoutImpl> base;

    std::vector<Point<int, GridLayoutImpl::dimension>> iCellForCentering;
    std::vector<std::array<double, GridLayoutImpl::dimension>> expectedPosition;

    std::vector<std::array<double, GridLayoutImpl::dimension>> actualPosition;


    void init()
    {
        auto& field           = base.field;
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        field = base.makeMyField_(layout->allocSize(currentQuantity));

        for (auto&& iCell : iCellForCentering)
        {
            auto const pos = base.layout->cellCenteredCoordinates(iCell);
            actualPosition.push_back(*pos);
        }
    }
};


template<typename GridLayoutImpl>
auto createCellCenteringParam()
{
    std::vector<GridLayoutCellCenteringParam<GridLayoutImpl>> params;

    std::string summaryName{"centeredCoords_summary"};
    std::string valueName{"centeredCoords_values"};

    std::string path{"./"};

    std::string summaryPath{path + summaryName + "_" + std::to_string(GridLayoutImpl::dimension)
                            + "d_O" + std::to_string(GridLayoutImpl::interp_order) + ".txt"};
    std::string valuePath{path + valueName + "_" + std::to_string(GridLayoutImpl::dimension) + "d_O"
                          + std::to_string(GridLayoutImpl::interp_order) + ".txt"};

    std::ifstream summary{summaryPath};
    std::ifstream value{valuePath};

    std::array<std::uint32_t, GridLayoutImpl::dimension> nbCell;
    std::array<double, GridLayoutImpl::dimension> dl;
    std::array<std::uint32_t, GridLayoutImpl::dimension> iStart;
    std::array<std::uint32_t, GridLayoutImpl::dimension> iEnd;
    std::array<double, GridLayoutImpl::dimension> origin;

    writeToArray(summary, nbCell);
    writeToArray(summary, dl);
    writeToArray(summary, iStart);
    writeToArray(summary, iEnd);
    writeToArray(summary, origin);

    params.emplace_back();

    params.back().base = createParam<GridLayoutImpl>(dl, nbCell, Point{origin});
    auto const& layout = *params.back().base.layout;


    std::array<std::uint32_t, GridLayoutImpl::dimension> icell;
    std::array<double, GridLayoutImpl::dimension> realPosition;

    while (writeToArray(value, icell) && writeToArray(value, realPosition))
    {
        auto& param = params[0];
        param.iCellForCentering.push_back(icell);
        param.expectedPosition.push_back(realPosition);
    }

    for (auto&& param : params)
    {
        param.init();
    }
    return params;
}



#endif
