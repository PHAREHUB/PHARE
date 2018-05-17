#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_FIELD_CENTERED_COORD_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_FIELD_CENTERED_COORD_H

#include <array>

#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "utilities/point/point.h"

namespace PHARE
{
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder>
struct GridLayoutFieldCenteringParam
{
    GridLayoutTestParam<GridLayoutImpl, dim, interpOrder> base;

    std::vector<std::array<uint32, dim>> iCellForCentering;
    std::vector<std::array<double, dim>> expectedPosition;

    std::vector<std::array<double, dim>> actualPosition;

    template<typename Array, std::size_t... I>
    auto fieldCoord_impl(Array const& array, std::index_sequence<I...>)
    {
        auto& field  = base.field;
        auto& layout = base.layout;
        auto& origin = base.origin;

        return layout->fieldNodeCoordinates(*field, origin, array[I]...);
    }

    template<typename T, std::size_t N, typename Indices = std::make_index_sequence<N>>
    auto fieldCoord(const std::array<T, N>& array)
    {
        return fieldCoord_impl(array, Indices{});
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
            pos = fieldCoord(iCell);

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
auto createFieldCenteringParam()
{
    std::vector<GridLayoutFieldCenteringParam<GridLayoutImpl, dim, interpOrder>> params;

    std::string summaryName{"fieldCoords_summary"};
    std::string valueName{"fieldCoords_values"};

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


    constexpr uint32 numberOfQuantities{14};

    while (!summary.eof())
    {
        std::string quantity;

        std::array<uint32, dim> nbCell;
        std::array<double, dim> dl;

        std::array<uint32, dim> iStart;
        std::array<uint32, dim> iEnd;

        std::array<double, dim> origin;

        summary >> quantity;

        if (summary.eof() || summary.bad())
            break;

        writeToArray(summary, nbCell);
        writeToArray(summary, dl);
        writeToArray(summary, iStart);
        writeToArray(summary, iEnd);
        writeToArray(summary, origin);

        params.emplace_back();

        // NOTE: before c++17 Point{origin} cannot deduce the corect type
        params.back().base
            = createParam<GridLayoutImpl, dim, interpOrder>(dl, nbCell, Point<double, dim>{origin});

        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
            params.back().base.currentQuantity = quantityIt->second;
    }



    while (!value.eof())
    {
        std::string quantity;
        std::array<uint32, dim> icell;
        std::array<double, dim> realPosition;

        value >> quantity;

        if (value.eof() || value.bad())
            break;

        writeToArray(value, icell);
        writeToArray(value, realPosition);

        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
        {
            auto hqIndex = static_cast<int>(quantityIt->second);

            auto& param = params[hqIndex];

            param.iCellForCentering.push_back(icell);
            param.expectedPosition.push_back(realPosition);
        }
    }

    for (auto&& param : params)
    {
        param.init();
    }
    return params;
}

} // namespace PHARE
#endif
