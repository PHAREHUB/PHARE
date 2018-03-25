#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <core/data/ndarray/ndarray_vector.h>
#include <string>


using PHARE::NdArrayVector1D;
using PHARE::NdArrayVector2D;
using PHARE::NdArrayVector3D;



TEST(NdArray1DTest, SizeIsOkAfterNonEmptyCreation)
{
    auto size = 10;
    NdArrayVector1D<> array1d{size};
    EXPECT_EQ(size, array1d.size());
}

TEST(NdArray2DTest, SizeIsOkAfterNonEmptyCreation)
{
    auto nx   = 12;
    auto ny   = 24;
    auto size = nx * ny;
    NdArrayVector2D<> array2d{nx, ny};
    EXPECT_EQ(size, array2d.size());
}


TEST(NdArray3DTest, SizeIsOkAfterNonEmptyCreation)
{
    auto nx   = 12;
    auto ny   = 24;
    auto nz   = 103;
    auto size = nx * ny * nz;
    NdArrayVector3D<> array3d{nx, ny, nz};
    EXPECT_EQ(size, array3d.size());
}


TEST(NdArray1DTest, IsModifiable)
{
    auto nx = 12;
    NdArrayVector1D<> array1d{nx};
    array1d(10) = 24;
    EXPECT_EQ(24, array1d(10));
}

TEST(NdArray1DTest, ArrayCanBeReadOnly)
{
    auto nx = 12;
    NdArrayVector1D<> array1d{nx};
    array1d(10)                  = 24;
    NdArrayVector1D<> const& ref = array1d;
    EXPECT_EQ(24, ref(10));
}


TEST(NdArray2DTest, IsModifiable)
{
    auto nx = 12;
    auto ny = 32;
    NdArrayVector2D<> array2d{nx, ny};
    array2d(10, 8) = 24;
    EXPECT_EQ(24, array2d(10, 8));
}


TEST(NdArray2DTest, ArrayCanBeReadOnly)
{
    auto nx = 12;
    auto ny = 87;
    NdArrayVector2D<> array2d{nx, ny};
    array2d(10, 2) = 24;
    NdArrayVector2D<> const& ref = array2d;
    EXPECT_EQ(24, ref(10, 2));
}


TEST(NdArray3DTest, IsModifiable)
{
    auto nx = 12;
    auto ny = 32;
    auto nz = 13;
    NdArrayVector3D<> array3d{nx, ny, nz};
    array3d(10, 8, 5) = 24;
    EXPECT_EQ(24, array3d(10, 8, 5));
}


TEST(NdArray3DTest, ArrayCanBeReadOnly)
{
    auto nx = 12;
    auto ny = 32;
    auto nz = 13;
    NdArrayVector3D<> array3d{nx, ny, nz};
    array3d(10, 8, 5) = 24;
    NdArrayVector3D<> const& ref = array3d;
    EXPECT_EQ(24, ref(10, 8, 5));
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
