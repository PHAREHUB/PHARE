#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_FIELD_CENTERED_COORD_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_FIELD_CENTERED_COORD_H

#include <array>

#include "core/data/grid/gridlayout.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/utilities/point/point.hpp"

using namespace PHARE::core;

template<typename GridLayoutImpl>
struct GridLayoutFieldCenteringParam
{
    GridLayoutTestParam<GridLayoutImpl> base;

    std::vector<std::array<std::uint32_t, GridLayoutImpl::dimension>> iCellForCentering;
    std::vector<std::array<double, GridLayoutImpl::dimension>> expectedPosition;
    std::vector<std::array<double, GridLayoutImpl::dimension>> actualPosition;

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
            Point<double, GridLayoutImpl::dimension> pos;
            pos = fieldCoord(iCell);

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
auto createFieldCenteringParam()
{
    std::vector<GridLayoutFieldCenteringParam<GridLayoutImpl>> params;

    std::string summaryName{"fieldCoords_summary"};
    std::string valueName{"fieldCoords_values"};

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


    // constexpr std::uint32_t numberOfQuantities{14};

    while (!summary.eof())
    {
        std::string quantity;

        std::array<std::uint32_t, GridLayoutImpl::dimension> nbCell;
        std::array<double, GridLayoutImpl::dimension> dl;
        std::array<std::uint32_t, GridLayoutImpl::dimension> iStart;
        std::array<std::uint32_t, GridLayoutImpl::dimension> iEnd;
        std::array<double, GridLayoutImpl::dimension> origin;

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
        params.back().base = createParam<GridLayoutImpl>(
            dl, nbCell, Point<double, GridLayoutImpl::dimension>{origin});

        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
            params.back().base.currentQuantity = quantityIt->second;
    }



    while (!value.eof())
    {
        std::string quantity;
        std::array<std::uint32_t, GridLayoutImpl::dimension> icell;
        std::array<double, GridLayoutImpl::dimension> realPosition;

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


#endif
