#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H


#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "hybrid/hybrid_quantities.h"
#include "utilities/box/box.h"
#include "utilities/point/point.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE;

std::vector<double> read(std::string filename);



template<typename GridLayoutImpl>
class a1DDerivative : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector1D<>, PHARE::HybridQuantity::Scalar> Ez;

public:
    a1DDerivative()
        : layout{{{0.1}}, {50}, Point{0.}, Box{Point{0}, Point{49}}}
        , By{"By", PHARE::HybridQuantity::Scalar::By,
             layout.allocSize(PHARE::HybridQuantity::Scalar::By)}
        , Ez{"Ez", PHARE::HybridQuantity::Scalar::Ez,
             layout.allocSize(PHARE::HybridQuantity::Scalar::Ez)}
    {
    }
};



template<typename GridLayoutImpl>
class a2DDerivative : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector2D<>, PHARE::HybridQuantity::Scalar> Ez;

public:
    a2DDerivative()
        : layout{{{0.1, 0.2}}, {50, 30}, Point{0., 0.}, Box{Point{0, 0}, Point{49, 39}}}
        , By{"By", PHARE::HybridQuantity::Scalar::By,
             layout.allocSize(PHARE::HybridQuantity::Scalar::By)}
        , Ez{"Ez", PHARE::HybridQuantity::Scalar::Ez,
             layout.allocSize(PHARE::HybridQuantity::Scalar::Ez)}
    {
    }
};



template<typename GridLayoutImpl>
class a3DDerivative : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> By;
    Field<NdArrayVector3D<>, PHARE::HybridQuantity::Scalar> Ez;

public:
    a3DDerivative()
        : layout{{{0.1, 0.2, 0.3}},
                 {50, 30, 40},
                 Point{0., 0., 0.},
                 Box{Point{0, 0, 0}, Point{49, 29, 39}}}
        , By{"By", PHARE::HybridQuantity::Scalar::By,
             layout.allocSize(PHARE::HybridQuantity::Scalar::By)}
        , Ez{"Ez", PHARE::HybridQuantity::Scalar::Ez,
             layout.allocSize(PHARE::HybridQuantity::Scalar::Ez)}
    {
    }
};


#endif
