#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_INDEXING_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_INDEXING_HPP

#include <array>
#include <fstream>

#include "core/data/grid/gridlayout.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/utilities/point/point.hpp"

using namespace PHARE::core;

template<auto options>
struct GridLayoutIndexingParam
{
    GridLayoutTestParam<options> base;

    std::array<std::uint32_t, options.dimension> actualPSI;
    std::array<std::uint32_t, options.dimension> actualPEI;
    std::array<std::uint32_t, options.dimension> actualGSI;
    std::array<std::uint32_t, options.dimension> actualGEI;

    std::array<std::uint32_t, options.dimension> expectedPSI;
    std::array<std::uint32_t, options.dimension> expectedPEI;
    std::array<std::uint32_t, options.dimension> expectedGSI;
    std::array<std::uint32_t, options.dimension> expectedGEI;


    void init()
    {
        auto& field           = base.field;
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        field = base.makeMyField_(layout->allocSize(currentQuantity));

        std::array<Direction, options.dimension> directions;
        directions[0] = Direction::X;
        if (options.dimension > 1)
        {
            directions[1] = Direction::Y;
        }
        if (options.dimension > 2)
        {
            directions[2] = Direction::Z;
        }

        for (std::uint32_t iDim = 0; iDim < options.dimension; ++iDim)
        {
            actualPSI[iDim] = layout->physicalStartIndex(*field, directions[iDim]);
            actualPEI[iDim] = layout->physicalEndIndex(*field, directions[iDim]);
            actualGSI[iDim] = layout->ghostStartIndex(*field, directions[iDim]);
            actualGEI[iDim] = layout->ghostEndIndex(*field, directions[iDim]);
        }
    }
};


template<auto options>
auto createIndexingParam()
{
    std::vector<GridLayoutIndexingParam<options>> params;

    std::string path{"./"};
    std::string baseName{"gridIndexing"};


    std::string fullName{path + baseName + "_" + std::to_string(options.dimension) + "d_O"
                         + std::to_string(options.interp_order) + ".txt"};
    std::ifstream inputFile{fullName};



    std::string layoutName{"yee"}; // hard coded for now


    if (!inputFile.is_open())
    {
        throw std::runtime_error("Error cannot open " + fullName);
    }

    std::uint32_t iQuantity;
    while (inputFile >> iQuantity)
    {
        std::array<std::uint32_t, options.dimension> numberCells;
        std::array<double, options.dimension> dl;

        writeToArray(inputFile, numberCells);
        writeToArray(inputFile, dl);

        params.emplace_back();
        params.back().base
            = createParam<options>(dl, numberCells, Point<double, options.dimension>{});

        writeToArray(inputFile, params.back().expectedPSI);
        writeToArray(inputFile, params.back().expectedPEI);
        writeToArray(inputFile, params.back().expectedGSI);
        writeToArray(inputFile, params.back().expectedGEI);

        params.back().base.currentQuantity = getQuantity(iQuantity);
    }

    return params;
}



#endif
