#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_H

#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayoutimplyee.h"
#include "core/utilities/types.h"


#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

template<typename GridLayoutImpl, template<typename> typename Param>
class GridLayoutTest : public ::testing::TestWithParam<Param<GridLayoutImpl>>
{
public:
    GridLayoutTest() = default;

    virtual void TearDown() override {}
    virtual void SetUp() override { param = this->GetParam(); }

    Param<GridLayoutImpl> param;
};



#endif
