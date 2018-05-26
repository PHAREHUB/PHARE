
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"

using namespace PHARE;


template<typename GridLayoutImpl>
class NextPrevTest : public ::testing::Test
{
protected:
    using layoutType = GridLayout<GridLayoutImpl>;
};

using layoutImpls
    = ::testing::Types<GridLayoutImplYee<1, 1>, GridLayoutImplYee<1, 2>, GridLayoutImplYee<1, 3>>;

TYPED_TEST_CASE(NextPrevTest, layoutImpls);




TYPED_TEST(NextPrevTest, nextPrimalIndex)
{
    auto nextPrimal = TestFixture::layoutType::nextIndex(QtyCentering::primal, 10);
    if (TestFixture::layoutType::nbrGhosts(QtyCentering::dual)
        > TestFixture::layoutType::nbrGhosts(QtyCentering::primal))
    {
        EXPECT_EQ(10, nextPrimal);
    }
    else
    {
        EXPECT_EQ(11, nextPrimal);
    }
}




TYPED_TEST(NextPrevTest, nextDualIndex)
{
    auto nextDual = TestFixture::layoutType::nextIndex(QtyCentering::dual, 10);
    if (TestFixture::layoutType::nbrGhosts(QtyCentering::dual)
        > TestFixture::layoutType::nbrGhosts(QtyCentering::primal))
    {
        EXPECT_EQ(11, nextDual);
    }
    else
    {
        EXPECT_EQ(10, nextDual);
    }
}



TYPED_TEST(NextPrevTest, prevPrimalIndex)
{
    auto prevDual = TestFixture::layoutType::prevIndex(QtyCentering::primal, 10);
    if (TestFixture::layoutType::nbrGhosts(QtyCentering::dual)
        > TestFixture::layoutType::nbrGhosts(QtyCentering::primal))
    {
        EXPECT_EQ(9, prevDual);
    }
    else
    {
        EXPECT_EQ(10, prevDual);
    }
}



TYPED_TEST(NextPrevTest, prevDualIndex)
{
    auto prevDual = TestFixture::layoutType::prevIndex(QtyCentering::dual, 10);
    if (TestFixture::layoutType::nbrGhosts(QtyCentering::dual)
        > TestFixture::layoutType::nbrGhosts(QtyCentering::primal))
    {
        EXPECT_EQ(10, prevDual);
    }
    else
    {
        EXPECT_EQ(9, prevDual);
    }
}
