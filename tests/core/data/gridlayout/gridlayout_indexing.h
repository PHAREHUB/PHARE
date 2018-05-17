#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_INDEXING_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_INDEXING_H

#include <array>
#include <fstream>

#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"

namespace PHARE
{
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
struct GridLayoutIndexingParam
{
    GridLayoutTestParam<GridLayoutImpl, dim, interpOrder> base;

    std::array<uint32, dim> actualPSI;
    std::array<uint32, dim> actualPEI;
    std::array<uint32, dim> actualGSI;
    std::array<uint32, dim> actualGEI;

    std::array<uint32, dim> expectedPSI;
    std::array<uint32, dim> expectedPEI;
    std::array<uint32, dim> expectedGSI;
    std::array<uint32, dim> expectedGEI;


    void init()
    {
        auto &field           = base.field;
        auto &layout          = base.layout;
        auto &currentQuantity = base.currentQuantity;

        field = base.makeMyField_(layout->allocSize(currentQuantity));

        std::array<Direction, dim> directions;
        directions[0] = Direction::X;
        if (dim > 1)
        {
            directions[1] = Direction::Y;
        }
        if (dim > 2)
        {
            directions[2] = Direction::Z;
        }

        for (uint32 iDim = 0; iDim < dim; ++iDim)
        {
            actualPSI[iDim] = layout->physicalStartIndex(*field, directions[iDim]);
            actualPEI[iDim] = layout->physicalEndIndex(*field, directions[iDim]);
            actualGSI[iDim] = layout->ghostStartIndex(*field, directions[iDim]);
            actualGEI[iDim] = layout->ghostEndIndex(*field, directions[iDim]);
        }
    }
};


template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
auto createIndexingParam()
{
    std::vector<GridLayoutIndexingParam<GridLayoutImpl, dim, interpOrder>> params;

    std::string path{"./"};
    std::string baseName{"gridIndexing"};


    std::string fullName{path + baseName + "_" + std::to_string(dim) + "d_O"
                         + std::to_string(interpOrder) + ".txt"};
    std::ifstream inputFile{fullName};



    std::string layoutName{"yee"}; // hard coded for now


    if (!inputFile.is_open())
    {
        throw std::runtime_error("Error cannot open " + fullName);
    }

    while (!inputFile.eof())
    {
        uint32 iQuantity;
        std::array<uint32, dim> numberCells;
        std::array<double, dim> dl;

        inputFile >> iQuantity;

        if (inputFile.eof() || inputFile.bad())
            break;

        writeToArray(inputFile, numberCells);
        writeToArray(inputFile, dl);

        params.emplace_back();
        params.back().base
            = createParam<GridLayoutImpl, dim, interpOrder>(dl, numberCells, Point<double, dim>{});

        writeToArray(inputFile, params.back().expectedPSI);
        writeToArray(inputFile, params.back().expectedPEI);
        writeToArray(inputFile, params.back().expectedGSI);
        writeToArray(inputFile, params.back().expectedGEI);

        params.back().base.currentQuantity = getQuantity(iQuantity);
    }

    return params;
}

} // namespace PHARE

#endif
