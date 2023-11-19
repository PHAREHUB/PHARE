#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/utilities/types.hpp"

#include "test_gridlayout.hpp"

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
