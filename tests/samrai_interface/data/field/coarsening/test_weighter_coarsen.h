#ifndef PHARE_TEST_WEIGHTER_COARSEN_H
#define PHARE_TEST_WEIGHTER_COARSEN_H

#include <numeric>

#include "data/field/coarsening/field_coarsen_index_weight.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::DoubleEq;
using testing::DoubleNear;
using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr_interface;

using GridYee1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using Field1D     = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

struct AWeighterData
{
    std::shared_ptr<CoarsenWeighter> weight;
};

struct AWeighter : public testing::TestWithParam<AWeighterData>
{
    void SetUp() override { param = GetParam(); }


    AWeighterData param;
};


TEST_P(AWeighter, hasSumOfWeightEqualToOne)
{
    auto weight = param.weight->weights();

    double totalWeight = std::accumulate(std::begin(weight), std::end(weight), 0.);


    EXPECT_THAT(totalWeight, DoubleEq(1.0));
}

AWeighterData createWeighter(std::size_t nbrPoints)
{
    AWeighterData weightData;
    weightData.weight = std::make_shared<CoarsenWeighter>(nbrPoints);
    return weightData;
}

INSTANTIATE_TEST_CASE_P(TestWithMultipleWeightPointsThat, AWeighter,
                        testing::ValuesIn({createWeighter(2), createWeighter(3), createWeighter(4),
                                           createWeighter(5), createWeighter(6), createWeighter(7),
                                           createWeighter(8), createWeighter(9), createWeighter(10),
                                           createWeighter(11)}));

#endif
