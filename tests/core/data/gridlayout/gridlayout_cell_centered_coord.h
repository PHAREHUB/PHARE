
#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_H

#include <array>

#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"

namespace PHARE
{
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
struct GridLayoutCellCenteringParam
{
    GridLayoutTestParam<GridLayoutImpl, dim, interpOrder> base;

    std::vector<std::array<uint32, dim>> iCellForCentering;
    std::vector<std::array<double, dim>> expectedPosition;

    std::vector<std::array<double, dim>> actualPosition;

    template<typename Array, std::size_t... I>
    auto cellCenteredCoord_impl(Array const& array, std::index_sequence<I...>)
    {
        return base.layout->cellCenteredCoordinates(array[I]...);
    }

    template<typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
    auto cellCenteredCoord(const std::array<T, N>& array)
    {
        return cellCenteredCoord_impl(array, Indices{});
    }

    void init()
    {
        auto& field           = base.field;
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        field = base.makeMyField_(layout->allocSize(currentQuantity));

        for (auto&& iCell : iCellForCentering)
        {
            Point<double, dim> pos;
            pos = cellCenteredCoord(iCell);

            std::array<double, dim> actualPos;

            for (std::size_t iDim = 0; iDim < dim; ++iDim)
            {
                actualPos[iDim] = pos[iDim];
            }


            actualPosition.push_back(actualPos);
        }
    }
};


template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
auto createCellCenteringParam()
{
    std::vector<GridLayoutCellCenteringParam<GridLayoutImpl, dim, interpOrder>> params;

    std::string summaryName{"centeredCoords_summary"};
    std::string valueName{"centeredCoords_values"};

    std::string path{"./"};

    std::string summaryPath{path + summaryName + "_" + std::to_string(dim) + "d_O"
                            + std::to_string(interpOrder) + ".txt"};
    std::string valuePath{path + valueName + "_" + std::to_string(dim) + "d_O"
                          + std::to_string(interpOrder) + ".txt"};

    std::ifstream summary{summaryPath};
    std::ifstream value{valuePath};

    std::string layoutName{"yee"};

    const std::map<std::string, HybridQuantity::Scalar> namesToQuantity{
        {"Bx", HybridQuantity::Scalar::Bx}, {"By", HybridQuantity::Scalar::By},
        {"Bz", HybridQuantity::Scalar::Bz}, {"Ex", HybridQuantity::Scalar::Ex},
        {"Ey", HybridQuantity::Scalar::Ey}, {"Ez", HybridQuantity::Scalar::Ez},
        {"Jx", HybridQuantity::Scalar::Jx}, {"Jy", HybridQuantity::Scalar::Jy},
        {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
        {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
        {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};

    std::array<uint32, dim> nbCell;
    std::array<double, dim> dl;

    std::array<uint32, dim> iStart;
    std::array<uint32, dim> iEnd;

    std::array<double, dim> origin;

    writeToArray(summary, nbCell);
    writeToArray(summary, dl);
    writeToArray(summary, iStart);
    writeToArray(summary, iEnd);
    writeToArray(summary, origin);



    params.emplace_back();

    // NOTE: c++17 : Point{origin}, C++14 : Point<double, dim>{origin}
    params.back().base
        = createParam<GridLayoutImpl, dim, interpOrder>(dl, nbCell, Point<double, dim>{origin});



    while (!value.eof())
    {
        std::array<uint32, dim> icell;
        std::array<double, dim> realPosition;

        if (value.eof() || value.bad())
            break;

        writeToArray(value, icell);
        writeToArray(value, realPosition);

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


} // namespace PHARE
#endif
