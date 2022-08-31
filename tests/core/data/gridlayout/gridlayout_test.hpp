#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"
#include "core/utilities/types.hpp"


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


template<typename GridLayout>
class TestGridLayout : public GridLayout
{ // to expose a default constructor
public:
    auto static constexpr dim = GridLayout::dimension;

    TestGridLayout() = default;

    TestGridLayout(std::uint32_t cells)
        : GridLayout{PHARE::core::ConstArray<double, dim>(1.0 / cells),
                     PHARE::core::ConstArray<std::uint32_t, dim>(cells),
                     PHARE::core::Point<double, dim>{PHARE::core::ConstArray<double, dim>(0)}}
    {
    }

    auto static make(std::uint32_t cells) { return TestGridLayout{cells}; }
};

#endif
