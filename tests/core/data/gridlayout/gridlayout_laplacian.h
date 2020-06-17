#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_LAPLACIAN_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_LAPLACIAN_H


#include "core/data/field/field.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "core/hybrid/hybrid_quantities.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE;

std::vector<double> read(std::string filename);



template<typename GridLayoutImpl>
class a1DLaplacian : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<1>, HybridQuantity::Scalar> Jz;

public:
    a1DLaplacian()
        : layout{{{0.1}}, {50}, Point<double, 1>{0.}}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
    {
    }
};



template<typename GridLayoutImpl>
class a2DLaplacian : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<2>, HybridQuantity::Scalar> Jz;

public:
    a2DLaplacian()
        : layout{{{0.1, 0.2}}, {50, 30}, Point<double, 2>{0., 0.}}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
    {
    }
};



template<typename GridLayoutImpl>
class a3DLaplacian : public ::testing::Test
{
protected:
    GridLayout<GridLayoutImpl> layout;
    static constexpr std::size_t interp_order = GridLayoutImpl::interp_order;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jx;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jy;
    Field<NdArrayVector<3>, HybridQuantity::Scalar> Jz;

public:
    a3DLaplacian()
        : layout{{{0.1, 0.2, 0.3}}, {50, 30, 40}, Point<double, 3>{0., 0., 0.}}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
    {
    }
};



#endif
