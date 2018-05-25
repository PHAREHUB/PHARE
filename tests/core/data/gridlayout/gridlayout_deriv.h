#ifndef PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H
#define PHARE_TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_DERIV_H

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/field/field.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "gridlayout_base_params.h"
#include "gridlayout_params.h"
#include "gridlayout_utilities.h"
#include "hybrid/hybrid_quantities.h"


using namespace PHARE;

std::vector<double> read(std::string filename)
{
    std::ifstream readFile(filename);
    assert(readFile.is_open());
    std::vector<double> x;

    std::copy(std::istream_iterator<double>(readFile), std::istream_iterator<double>(),
              std::back_inserter(x));
    return x;
}



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
        : layout{{{0.1}}, {50}, Point<double, 1>{0.}}
        , By{"By", PHARE::HybridQuantity::Scalar::By,
             layout.allocSize(PHARE::HybridQuantity::Scalar::By)}
        , Ez{"Ez", PHARE::HybridQuantity::Scalar::Ez,
             layout.allocSize(PHARE::HybridQuantity::Scalar::Ez)}
    {
    }
};




#endif
