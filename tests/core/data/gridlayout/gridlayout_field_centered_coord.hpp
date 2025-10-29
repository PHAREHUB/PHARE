#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_FIELD_CENTERED_COORD_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_FIELD_CENTERED_COORD_HPP


#include "core/utilities/point/point.hpp"

#include "gridlayout_params.hpp"
#include "gridlayout_base_params.hpp"

#include <array>

using namespace PHARE::core;

template<typename GridLayoutImpl>
struct GridLayoutFieldCenteringParam
{
    GridLayoutTestParam<GridLayoutImpl> base;

    std::vector<Point<int, GridLayoutImpl::dimension>> iCellForCentering;
    std::vector<std::array<double, GridLayoutImpl::dimension>> expectedPosition;
    std::vector<std::array<double, GridLayoutImpl::dimension>> actualPosition;


    void init()
    {
        auto& field           = base.field;
        auto& layout          = base.layout;
        auto& currentQuantity = base.currentQuantity;

        field = base.makeMyField_(layout->allocSize(currentQuantity));

        for (auto&& iCell : iCellForCentering)
        {
            auto const pos = layout->fieldNodeCoordinates(*field, iCell);
            actualPosition.push_back(*pos);
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

    std::map<std::string, HybridQuantity::Scalar> const namesToQuantity{
        {"Bx", HybridQuantity::Scalar::Bx}, {"By", HybridQuantity::Scalar::By},
        {"Bz", HybridQuantity::Scalar::Bz}, {"Ex", HybridQuantity::Scalar::Ex},
        {"Ey", HybridQuantity::Scalar::Ey}, {"Ez", HybridQuantity::Scalar::Ez},
        {"Jx", HybridQuantity::Scalar::Jx}, {"Jy", HybridQuantity::Scalar::Jy},
        {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
        {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
        {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};


    // constexpr std::uint32_t numberOfQuantities{14};

    std::string quantity;
    while (summary >> quantity)
    {
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


        params.back().base = createParam<GridLayoutImpl>(dl, nbCell, Point{origin});

        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
            params.back().base.currentQuantity = quantityIt->second;
    }


    std::array<std::uint32_t, GridLayoutImpl::dimension> icell;
    std::array<double, GridLayoutImpl::dimension> realPosition;

    while (value >> quantity && writeToArray(value, icell) && writeToArray(value, realPosition))
    {
        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
        {
            auto hqIndex = static_cast<std::size_t>(quantityIt->second);

            assert(hqIndex < params.size());
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
