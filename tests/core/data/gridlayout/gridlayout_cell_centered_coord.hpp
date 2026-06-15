
#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_HPP

#include <array>

#include "core/data/grid/gridlayout.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/utilities/point/point.hpp"

using namespace PHARE::core;


template<auto options>
struct GridLayoutCellCenteringParam
{
    GridLayoutTestParam<options> base;

    std::vector<Point<int, options.dimension>> iCellForCentering;
    std::vector<std::array<double, options.dimension>> expectedPosition;

    std::vector<std::array<double, options.dimension>> actualPosition;


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


template<auto options>
auto createCellCenteringParam()
{
    std::vector<GridLayoutCellCenteringParam<options>> params;

    std::string summaryName{"centeredCoords_summary"};
    std::string valueName{"centeredCoords_values"};

    std::string path{"./"};

    std::string summaryPath{path + summaryName + "_" + std::to_string(options.dimension) + "d_O"
                            + std::to_string(options.interp_order) + ".txt"};
    std::string valuePath{path + valueName + "_" + std::to_string(options.dimension) + "d_O"
                          + std::to_string(options.interp_order) + ".txt"};

    std::ifstream summary{summaryPath};
    std::ifstream value{valuePath};

    std::array<std::uint32_t, options.dimension> nbCell;
    std::array<double, options.dimension> dl;
    std::array<std::uint32_t, options.dimension> iStart;
    std::array<std::uint32_t, options.dimension> iEnd;
    std::array<double, options.dimension> origin;

    writeToArray(summary, nbCell);
    writeToArray(summary, dl);
    writeToArray(summary, iStart);
    writeToArray(summary, iEnd);
    writeToArray(summary, origin);

    params.emplace_back();

    params.back().base = createParam<options>(dl, nbCell, Point{origin});

    std::array<std::uint32_t, options.dimension> icell;
    std::array<double, options.dimension> realPosition;

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
