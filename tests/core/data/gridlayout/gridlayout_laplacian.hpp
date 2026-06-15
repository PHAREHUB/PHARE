#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_LAPLACIAN_HPP
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_LAPLACIAN_HPP


#include "core/data/grid/grid.hpp"
#include "core/utilities/point/point.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;

std::vector<double> read(std::string filename);



template<typename TestParam_t>
class a1DLaplacian : public ::testing::Test
{
    auto constexpr static options = TestParam_t::field_options;

protected:
    TestParam_t::GridLayout_t layout;
    static constexpr std::size_t interp_order = options.interp_order;
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> Jx;
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> Jy;
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> Jz;

public:
    a1DLaplacian()
        : layout{{{0.1}}, {50}, Point<double, 1>{0.}}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
    {
    }
};



template<typename TestParam_t>
class a2DLaplacian : public ::testing::Test
{
    auto constexpr static options = TestParam_t::field_options;

protected:
    TestParam_t::GridLayout_t layout;
    static constexpr std::size_t interp_order = options.interp_order;
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> Jx;
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> Jy;
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> Jz;

public:
    a2DLaplacian()
        : layout{{{0.1, 0.2}}, {50, 30}, Point<double, 2>{0., 0.}}
        , Jx{"Jx", HybridQuantity::Scalar::Jx, layout.allocSize(HybridQuantity::Scalar::Jx)}
        , Jy{"Jy", HybridQuantity::Scalar::Jy, layout.allocSize(HybridQuantity::Scalar::Jy)}
        , Jz{"Jz", HybridQuantity::Scalar::Jz, layout.allocSize(HybridQuantity::Scalar::Jz)}
    {
    }
};



template<typename TestParam_t>
class a3DLaplacian : public ::testing::Test
{
    auto constexpr static options = TestParam_t::field_options;

protected:
    TestParam_t::GridLayout_t layout;
    static constexpr std::size_t interp_order = options.interp_order;
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> Jx;
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> Jy;
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> Jz;

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
