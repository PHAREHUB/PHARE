#ifndef TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_H
#define TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_H


#include "core/data/grid/gridlayout.h"
#include "core/data/ndarray/ndarray_vector.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "core/utilities/point/point.h"

#include <fstream>
#include <vector>



using namespace PHARE::core;


template<typename GridLayoutImpl>
struct GridLayoutAllocSizeParam
{
    GridLayoutTestParam<GridLayoutImpl> base;

    std::array<uint32, GridLayoutImpl::dimension> expectedAllocSize;
    std::array<uint32, GridLayoutImpl::dimension> expectedAllocSizeDerived;

    std::array<uint32, GridLayoutImpl::dimension> actualAllocSize;
    std::array<uint32, GridLayoutImpl::dimension> actualAllocSizeDerived;


    void init()
    {
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        actualAllocSize = layout->allocSize(currentQuantity);

        actualAllocSizeDerived[0] = layout->allocSizeDerived(currentQuantity, Direction::X)[0];
        if (GridLayoutImpl::dimension > 1)
        {
            actualAllocSizeDerived[1] = layout->allocSizeDerived(currentQuantity, Direction::Y)[1];
        }
        if (GridLayoutImpl::dimension > 2)
        {
            actualAllocSizeDerived[2] = layout->allocSizeDerived(currentQuantity, Direction::Z)[2];
        }
    }
};


template<typename GridLayoutImpl>
auto createAllocSizeParam()
{
    std::vector<GridLayoutAllocSizeParam<GridLayoutImpl>> params;

    std::string path{"./"};
    std::string baseName{"allocSizes"};


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

        writeToArray(inputFile, params.back().expectedAllocSize);
        writeToArray(inputFile, params.back().expectedAllocSizeDerived);

        params.back().base.currentQuantity = getQuantity(iQuantity);
    }

    return params;
}



#endif
