#ifndef TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_H
#define TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_H


#include "data/grid/gridlayout.h"
#include "data/ndarray/ndarray_vector.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"

#include <fstream>
#include <vector>




namespace PHARE
{
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
struct GridLayoutAllocSizeParam
{
    GridLayoutTestParam<GridLayoutImpl, dim, interpOrder> base;

    std::array<uint32, dim> expectedAllocSize;
    std::array<uint32, dim> expectedAllocSizeDerived;

    std::array<uint32, dim> actualAllocSize;
    std::array<uint32, dim> actualAllocSizeDerived;


    void init()
    {
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        actualAllocSize = layout->allocSize(currentQuantity);

        actualAllocSizeDerived[0] = layout->allocSizeDerived(currentQuantity, Direction::X)[0];
        if (dim > 1)
        {
            actualAllocSizeDerived[1] = layout->allocSizeDerived(currentQuantity, Direction::Y)[1];
        }
        if (dim > 2)
        {
            actualAllocSizeDerived[2] = layout->allocSizeDerived(currentQuantity, Direction::Z)[2];
        }
    }
};


template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
auto createAllocSizeParam()
{
    std::vector<GridLayoutAllocSizeParam<GridLayoutImpl, dim, interpOrder>> params;

    std::string path{"./"};
    std::string baseName{"allocSizes"};


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

        writeToArray(inputFile, params.back().expectedAllocSize);
        writeToArray(inputFile, params.back().expectedAllocSizeDerived);

        params.back().base.currentQuantity = getQuantity(iQuantity);
    }

    return params;
}

} // namespace PHARE

#endif
