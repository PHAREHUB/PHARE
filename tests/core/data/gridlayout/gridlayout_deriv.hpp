#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H


#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "gridlayout_base_params.hpp"
#include "gridlayout_params.hpp"
#include "gridlayout_utilities.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

std::vector<double> read(std::string filename);



template<typename GridLayoutImpl>
class a1DDerivative : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Ez;

public:
    a1DDerivative()
        : layout{{{0.1}}, {50}, Point{0.}}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Ez{"Ez", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
    {
    }
};



template<typename GridLayoutImpl>
class a2DDerivative : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Ez;

public:
    a2DDerivative()
        : layout{{{0.1, 0.2}}, {50, 30}, Point{0., 0.}}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Ez{"Ez", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
    {
    }
};



template<typename GridLayoutImpl>
class a3DDerivative : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> By;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Ez;

public:
    a3DDerivative()
        : layout{{{0.1, 0.2, 0.3}}, {50, 30, 40}, Point{0., 0., 0.}}
        , By{"By", HybridQuantity::Scalar::By, layout.allocSize(HybridQuantity::Scalar::By)}
        , Ez{"Ez", HybridQuantity::Scalar::Ez, layout.allocSize(HybridQuantity::Scalar::Ez)}
    {
    }
};


#endif
