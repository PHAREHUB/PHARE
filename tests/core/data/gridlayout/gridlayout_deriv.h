#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H

#include "data/grid/gridlayout.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"

namespace PHARE
{
template<Layout layoutType, std::size_t dim>
struct GridLayoutDerivParam
{
    GridLayoutTestParam<layoutType, dim> base;

    std::vector<std::array<uint32, dim>> iCell;

    std::vector<std::array<uint32, dim>> iCellDeriv;

    std::vector<double> fieldValues;

    std::vector<double> expectedDeriv;

    decltype(base.field) derivedField;

    HybridQuantity::Scalar derivedQuantity;


    void init()
    {
        // only 1D for now
        auto &field           = base.field;
        auto &layout          = base.layout;
        auto &currentQuantity = base.currentQuantity;

        field = std::make_shared<Field<NdArrayVector1D<>, HybridQuantity::Scalar>>(
            "field", currentQuantity, layout->allocSize(currentQuantity)[0]);
        derivedField = std::make_shared<Field<NdArrayVector1D<>, HybridQuantity::Scalar>>(
            "derived", derivedQuantity, layout->allocSizeDerived(derivedQuantity, Direction::X)[0]);

        for (std::size_t index = 0; index < iCell.size(); ++index)
        {
            (*field)(iCell[index][0]) = fieldValues[index];
        }

        layout->deriv(*field, Direction::X, *derivedField);
    }
};



template<Layout layout, std::size_t dim>
auto createDerivParam()
{
    std::vector<GridLayoutDerivParam<layout, dim>> params;

    std::string summaryName{"deriv_summary"};
    std::string valueName{"deriv_values"};
    std::string derivedValue{"deriv_values_derived"};


    std::string path{"./"};

    std::string strDim = std::to_string(dim);

    std::string summaryPath{path + summaryName + "_" + strDim + "d.txt"};
    std::string valuePath{path + valueName + "_" + strDim + "d.txt"};
    std::string derivedPath{path + derivedValue + "_" + strDim + "d.txt"};

    std::ifstream summary{summaryPath};
    std::ifstream value{valuePath};
    std::ifstream derived{derivedPath};

    std::string layoutName{"yee"};

    const std::map<std::string, HybridQuantity::Scalar> namesToQuantity{
        {"Bx", HybridQuantity::Scalar::Bx}, {"By", HybridQuantity::Scalar::By},
        {"Bz", HybridQuantity::Scalar::Bz}, {"Ex", HybridQuantity::Scalar::Ex},
        {"Ey", HybridQuantity::Scalar::Ey}, {"Ez", HybridQuantity::Scalar::Ez},
        {"Jx", HybridQuantity::Scalar::Jx}, {"Jy", HybridQuantity::Scalar::Jy},
        {"Jz", HybridQuantity::Scalar::Jz}, {"rho", HybridQuantity::Scalar::rho},
        {"Vx", HybridQuantity::Scalar::Vx}, {"Vy", HybridQuantity::Scalar::Vy},
        {"Vz", HybridQuantity::Scalar::Vz}, {"P", HybridQuantity::Scalar::P}};

    const std::map<HybridQuantity::Scalar, std::size_t> quantityToIndex{
        {HybridQuantity::Scalar::Ex, 0},
        {HybridQuantity::Scalar::Ey, 1},
        {HybridQuantity::Scalar::rho, 0},
        {HybridQuantity::Scalar::Bz, 1}};


    constexpr uint32 numberOfQuantities{2};

    /* std::map<std::string, uint32> functionNameToIndex; */
    /* bool functionNameMapIsComplete{false}; */
    /* uint32 countFunctionName{0}; */

    while (!summary.eof())
    {
        int currentOrder{0};
        std::string quantity;
        std::string derivedQuantity;


        std::array<uint32, dim> nbCell;
        std::array<double, dim> dl;

        std::array<uint32, dim> iGhostStart;
        std::array<uint32, dim> iGhostEnd;

        std::array<uint32, dim> iStart;
        std::array<uint32, dim> iEnd;

        std::array<double, dim> origin;

        summary >> currentOrder;
        summary >> quantity;
        summary >> derivedQuantity;

        /* if (!functionNameMapIsComplete) */
        /* { */
        /*     if (functionNameToIndex.find(functionName) == functionNameToIndex.end()) */
        /*     { */
        /*         functionNameToIndex[functionName] = countFunctionName; */
        /*         ++countFunctionName; */
        /*     } */
        /*     else */
        /*     { */
        /*         functionNameMapIsComplete = true; */
        /*     } */
        /* } */

        if (summary.eof() || summary.bad())
            break;

        writeToArray(summary, nbCell);
        writeToArray(summary, dl);
        writeToArray(summary, iGhostStart);
        writeToArray(summary, iGhostEnd);
        writeToArray(summary, iStart);
        writeToArray(summary, iEnd);
        writeToArray(summary, origin);

        params.emplace_back();

        // NOTE: before c++17 Point{origin} cannot deduce the corect type
        params.back().base = createParam<layout, dim>(layoutName, currentOrder, dl, nbCell,
                                                      Point<double, dim>{origin});


        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
            params.back().base.currentQuantity = quantityIt->second;

        auto derivQuantityIt = namesToQuantity.find(derivedQuantity);
        if (derivQuantityIt != namesToQuantity.end())
            params.back().derivedQuantity = quantityIt->second;
    }



    while (!value.eof())
    {
        int order{0};
        std::string quantity;

        std::array<uint32, dim> icell;

        double functionValue;

        value >> order;
        value >> quantity;

        if (value.eof() || value.bad())
            break;

        writeToArray(value, icell);

        value >> functionValue;

        auto quantityIt = namesToQuantity.find(quantity);
        if (quantityIt != namesToQuantity.end())
        {
            std::size_t hqIndex = 0;
            auto hqIndexIt      = quantityToIndex.find(quantityIt->second);

            if (hqIndexIt != quantityToIndex.end())
                hqIndex = hqIndexIt->second;

            auto &param = params[(order - 1) * numberOfQuantities + hqIndex];


            param.iCell.push_back(icell);
            param.fieldValues.push_back(functionValue);
        }
    }
    while (!derived.eof())
    {
        int order{0};
        std::string derivQuantity;

        std::array<uint32, dim> icell;

        double functionDeriv;

        derived >> order;
        derived >> derivQuantity;

        if (value.eof() || value.bad())
            break;

        writeToArray(derived, icell);

        derived >> functionDeriv;

        auto derivQuantityIt = namesToQuantity.find(derivQuantity);
        if (derivQuantityIt != namesToQuantity.end())
        {
            std::size_t hqIndex = 0;
            auto hqIndexIt      = quantityToIndex.find(derivQuantityIt->second);

            if (hqIndexIt != quantityToIndex.end())
                hqIndex = hqIndexIt->second;

            /* auto functionIndexIt = functionNameToIndex.find(functionName); */

            /* if (functionIndexIt != functionNameToIndex.end()) */
            /* { */
            /*     std::size_t numberOfFunction = functionNameToIndex.size(); */
            auto &param = params[(order - 1) * numberOfQuantities + hqIndex];


            param.iCellDeriv.push_back(icell);
            param.expectedDeriv.push_back(functionDeriv);
            /* } */
        }
    }

    for (auto &&param : params)
    {
        param.init();
    }
    return params;
}



} // namespace PHARE

#endif
