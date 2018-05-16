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
template<typename GridLayoutImpl, std::size_t dim, std::size_t interpOrder,
         template<typename, std::size_t, std::size_t> typename Param>
class GridLayoutTest : public ::testing::TestWithParam<Param<GridLayoutImpl, dim, interpOrder>>
{
public:
    GridLayoutTest() = default;

    virtual void TearDown() override {}
    virtual void SetUp() override { param = this->GetParam(); }

    Param<GridLayoutImpl, dim, interpOrder> param;
};


} // namespace PHARE

#endif
