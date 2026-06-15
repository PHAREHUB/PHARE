#ifndef TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_HPP


#include "core/data/grid/gridlayout.hpp"
#include "core/utilities/point/point.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "gridlayout_base_params.hpp"

#include <fstream>
#include <vector>



using namespace PHARE::core;


template<auto options>
struct GridLayoutAllocSizeParam
{
    GridLayoutTestParam<options> base;

    std::array<std::uint32_t, options.dimension> expectedAllocSize;
    std::array<std::uint32_t, options.dimension> expectedAllocSizeDerived;

    std::array<std::uint32_t, options.dimension> actualAllocSize;
    std::array<std::uint32_t, options.dimension> actualAllocSizeDerived;


    void init()
    {
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        actualAllocSize = layout->allocSize(currentQuantity);

        actualAllocSizeDerived[0] = layout->allocSizeDerived(currentQuantity, Direction::X)[0];
        if (options.dimension > 1)
        {
            actualAllocSizeDerived[1] = layout->allocSizeDerived(currentQuantity, Direction::Y)[1];
        }
        if (options.dimension > 2)
        {
            actualAllocSizeDerived[2] = layout->allocSizeDerived(currentQuantity, Direction::Z)[2];
        }
    }
};


template<auto options>
auto createAllocSizeParam()
{
    std::vector<GridLayoutAllocSizeParam<options>> params;

    std::string path{"./"};
    std::string baseName{"allocSizes"};


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

        writeToArray(inputFile, params.back().expectedAllocSize);
        writeToArray(inputFile, params.back().expectedAllocSizeDerived);

        params.back().base.currentQuantity = getQuantity(iQuantity);
    }

    return params;
}



#endif
