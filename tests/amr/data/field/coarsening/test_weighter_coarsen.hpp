#ifndef PHARE_TEST_WEIGHTER_COARSEN_HPP
#define PHARE_TEST_WEIGHTER_COARSEN_HPP


#include "amr/data/field/coarsening/coarsen_weighter.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <numeric>

using testing::DoubleEq;
using testing::DoubleNear;
using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;


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
    auto weight        = param.weight->weights();
    double totalWeight = std::accumulate(std::begin(weight), std::end(weight), 0.);
    EXPECT_THAT(totalWeight, DoubleEq(1.0));
}


AWeighterData createWeighter(std::size_t nbrPoints)
{
    AWeighterData weightData;
    weightData.weight = std::make_shared<CoarsenWeighter>(nbrPoints);
    return weightData;
}


INSTANTIATE_TEST_SUITE_P(TestWithMultipleWeightPointsThat, AWeighter,
                         testing::ValuesIn({createWeighter(2), createWeighter(3), createWeighter(4),
                                            createWeighter(5), createWeighter(6), createWeighter(7),
                                            createWeighter(8), createWeighter(9),
                                            createWeighter(10), createWeighter(11)}));

#endif
