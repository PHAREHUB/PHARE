#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_H
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_H

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayoutimplyee.h"
#include "utilities/types.h"

/* #include "gridlayout_params.h" */

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace PHARE
{
template<typename GridLayoutImpl, template<typename> typename Param>
class GridLayoutTest : public ::testing::TestWithParam<Param<GridLayoutImpl>>
{
public:
    GridLayoutTest() = default;

    virtual void TearDown() override {}
    virtual void SetUp() override { param = this->GetParam(); }

    Param<GridLayoutImpl> param;
};


} // namespace PHARE

#endif
