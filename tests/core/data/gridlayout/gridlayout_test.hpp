#ifndef TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP
#define TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP


#include "gtest/gtest.h"


template<auto options, template<auto> typename Param>
class GridLayoutTest : public ::testing::TestWithParam<Param<options>>
{
public:
    using Options = decltype(options);

    GridLayoutTest() = default;

    virtual void TearDown() override {}
    virtual void SetUp() override { param = this->GetParam(); }

    Param<options> param;
};


#endif // TESTS_CORE_DATA_GRIDLAYOUT_GRIDLAYOUT_TEST_HPP
