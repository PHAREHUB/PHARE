#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <memory>

#include "core/numerics/ampere/ampere.h"

using namespace PHARE;

struct GridLayoutMock1D
{
    static const int dimension = 1;
    void deriv() {}
};

struct GridLayoutMock2D
{
    static const int dimension = 2;
};

struct GridLayoutMock3D
{
    static const int dimension = 3;
};

struct VecFieldMock
{
};

TEST(Ampere, canBe1D)
{
    Ampere<GridLayoutMock1D> ampere;
}


TEST(Ampere, shouldBeGivenAGridLayoutPointerToBeOperational)
{
    Ampere<GridLayoutMock1D> ampere;
    auto layout = std::make_unique<GridLayoutMock1D>();
    VecFieldMock B, J;
    EXPECT_ANY_THROW(ampere(B, J));
    ampere.setLayout(layout.get());
    EXPECT_NO_THROW(ampere(B, J));
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
