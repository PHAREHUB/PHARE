
#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_CELL_CENTERED_COORD_H

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

    std::vector<std::array<std::uint32_t, GridLayoutImpl::dimension>> iCellForCentering;
    std::vector<std::array<double, GridLayoutImpl::dimension>> expectedPosition;

    std::vector<std::array<double, GridLayoutImpl::dimension>> actualPosition;

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
            Point<double, GridLayoutImpl::dimension> pos;
            pos = cellCenteredCoord(iCell);

            std::array<double, GridLayoutImpl::dimension> actualPos;

            for (std::size_t iDim = 0; iDim < GridLayoutImpl::dimension; ++iDim)
            {
                actualPos[iDim] = pos[iDim];
            }


            actualPosition.push_back(actualPos);
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

    std::string layoutName{"yee"};

    const std::map<std::string, HybridQuantity::Scalar> namesToQuantity{
        {"Bx", HybridQuantity::Scalar::Bx}, {"By", HybridQuantity::Scalar::By},
        {"Bz", HybridQuantity::Scalar::Bz}, {"Ex", HybridQuantity::Scalar::Ex},
        {"Ey", HybridQuantity::Scalar::Ey}, {"Ez", HybridQuantity::Scalar::Ez},
        {"Jx", HybridQuantity::Scalar::Jx}, {"Jy", HybridQuantity::Scalar::Jy},
        {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
        {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
        {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};

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

    // NOTE: c++17 : Point{origin}, C++14 : Point<double, dim>{origin}
    params.back().base
        = createParam<GridLayoutImpl>(dl, nbCell, Point<double, GridLayoutImpl::dimension>{origin});



    while (!value.eof())
    {
        std::array<std::uint32_t, GridLayoutImpl::dimension> icell;
        std::array<double, GridLayoutImpl::dimension> realPosition;

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



#endif
