#ifndef TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_ALLOCSIZE_HPP


#include "core/data/grid/gridlayout.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/utilities/point/point.hpp"

#include <fstream>
#include <vector>



using namespace PHARE::core;


template<typename GridLayoutImpl>
struct GridLayoutAllocSizeParam
{
    GridLayoutTestParam<GridLayoutImpl> base;

    std::array<std::uint32_t, GridLayoutImpl::dimension> expectedAllocSize;
    std::array<std::uint32_t, GridLayoutImpl::dimension> expectedAllocSizeDerived;

    std::array<std::uint32_t, GridLayoutImpl::dimension> actualAllocSize;
    std::array<std::uint32_t, GridLayoutImpl::dimension> actualAllocSizeDerived;


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

    std::uint32_t iQuantity;
    while (inputFile >> iQuantity)
    {
        std::array<std::uint32_t, GridLayoutImpl::dimension> numberCells;
        std::array<double, GridLayoutImpl::dimension> dl;

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
