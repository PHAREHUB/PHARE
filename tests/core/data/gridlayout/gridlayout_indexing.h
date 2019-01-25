#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_INDEXING_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_INDEXING_H

#include <array>
#include <fstream>

#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"

using namespace PHARE::core;

template<typename GridLayoutImpl>
struct GridLayoutIndexingParam
{
    GridLayoutTestParam<GridLayoutImpl> base;

    std::array<uint32, GridLayoutImpl::dimension> actualPSI;
    std::array<uint32, GridLayoutImpl::dimension> actualPEI;
    std::array<uint32, GridLayoutImpl::dimension> actualGSI;
    std::array<uint32, GridLayoutImpl::dimension> actualGEI;

    std::array<uint32, GridLayoutImpl::dimension> expectedPSI;
    std::array<uint32, GridLayoutImpl::dimension> expectedPEI;
    std::array<uint32, GridLayoutImpl::dimension> expectedGSI;
    std::array<uint32, GridLayoutImpl::dimension> expectedGEI;


    void init()
    {
        auto& field           = base.field;
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        field = base.makeMyField_(layout->allocSize(currentQuantity));

        std::array<Direction, GridLayoutImpl::dimension> directions;
        directions[0] = Direction::X;
        if (GridLayoutImpl::dimension > 1)
        {
            directions[1] = Direction::Y;
        }
        if (GridLayoutImpl::dimension > 2)
        {
            directions[2] = Direction::Z;
        }

        for (uint32 iDim = 0; iDim < GridLayoutImpl::dimension; ++iDim)
        {
            actualPSI[iDim] = layout->physicalStartIndex(*field, directions[iDim]);
            actualPEI[iDim] = layout->physicalEndIndex(*field, directions[iDim]);
            actualGSI[iDim] = layout->ghostStartIndex(*field, directions[iDim]);
            actualGEI[iDim] = layout->ghostEndIndex(*field, directions[iDim]);
        }
    }
};


template<typename GridLayoutImpl>
auto createIndexingParam()
{
    std::vector<GridLayoutIndexingParam<GridLayoutImpl>> params;

    std::string path{"./"};
    std::string baseName{"gridIndexing"};


    std::string fullName{path + baseName + "_" + std::to_string(GridLayoutImpl::dimension) + "d_O"
                         + std::to_string(GridLayoutImpl::interp_order) + ".txt"};
    std::ifstream inputFile{fullName};



    std::string layoutName{"yee"}; // hard coded for now


    if (!inputFile.is_open())
    {
        throw std::runtime_error("Error cannot open " + fullName);
    }

    while (!inputFile.eof())
    {
        uint32 iQuantity;
        std::array<uint32, GridLayoutImpl::dimension> numberCells;
        std::array<double, GridLayoutImpl::dimension> dl;

        inputFile >> iQuantity;

        if (inputFile.eof() || inputFile.bad())
            break;

        writeToArray(inputFile, numberCells);
        writeToArray(inputFile, dl);

        params.emplace_back();
        params.back().base = createParam<GridLayoutImpl>(
            dl, numberCells, Point<double, GridLayoutImpl::dimension>{});

        writeToArray(inputFile, params.back().expectedPSI);
        writeToArray(inputFile, params.back().expectedPEI);
        writeToArray(inputFile, params.back().expectedGSI);
        writeToArray(inputFile, params.back().expectedGEI);

        params.back().base.currentQuantity = getQuantity(iQuantity);
    }

    return params;
}



#endif
