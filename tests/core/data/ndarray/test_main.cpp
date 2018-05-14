#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <random>
#include <string>

#include <core/data/ndarray/ndarray_vector.h>


using PHARE::NdArrayVector1D;
using PHARE::NdArrayVector2D;
using PHARE::NdArrayVector3D;



TEST(NdArray1DTest, SizeIsOkAfterNonEmptyCreation)
{
    auto size = 10u;
    NdArrayVector1D<> array1d{size};
    EXPECT_EQ(size, array1d.size());
}

TEST(NdArray2DTest, SizeIsOkAfterNonEmptyCreation)
{
    auto nx   = 12u;
    auto ny   = 24u;
    auto size = nx * ny;
    NdArrayVector2D<> array2d{nx, ny};
    EXPECT_EQ(size, array2d.size());
}


TEST(NdArray3DTest, SizeIsOkAfterNonEmptyCreation)
{
    auto nx   = 12u;
    auto ny   = 24u;
    auto nz   = 103u;
    auto size = nx * ny * nz;
    NdArrayVector3D<> array3d{nx, ny, nz};
    EXPECT_EQ(size, array3d.size());
}


TEST(NdArray1DTest, IsModifiable)
{
    auto nx = 12u;
    NdArrayVector1D<> array1d{nx};
    array1d(10) = 24;
    EXPECT_EQ(24, array1d(10));
}

TEST(NdArray1DTest, ArrayCanBeReadOnly)
{
    auto nx = 12u;
    NdArrayVector1D<> array1d{nx};
    array1d(10)                  = 24;
    NdArrayVector1D<> const& ref = array1d;
    EXPECT_EQ(24, ref(10));
}


TEST(NdArray2DTest, IsModifiable)
{
    auto nx = 12u;
    auto ny = 32u;
    NdArrayVector2D<> array2d{nx, ny};
    array2d(10, 8) = 24;
    EXPECT_EQ(24, array2d(10, 8));
}


TEST(NdArray2DTest, ArrayCanBeReadOnly)
{
    auto nx = 12u;
    auto ny = 87u;
    NdArrayVector2D<> array2d{nx, ny};
    array2d(10, 2)               = 24;
    NdArrayVector2D<> const& ref = array2d;
    EXPECT_EQ(24, ref(10, 2));
}


TEST(NdArray3DTest, IsModifiable)
{
    auto nx = 12u;
    auto ny = 32u;
    auto nz = 13u;
    NdArrayVector3D<> array3d{nx, ny, nz};
    array3d(10, 8, 5) = 24;
    EXPECT_EQ(24, array3d(10, 8, 5));
}


TEST(NdArray3DTest, ArrayCanBeReadOnly)
{
    auto nx = 12u;
    auto ny = 32u;
    auto nz = 13u;
    NdArrayVector3D<> array3d{nx, ny, nz};
    array3d(10, 8, 5)            = 24;
    NdArrayVector3D<> const& ref = array3d;
    EXPECT_EQ(24, ref(10, 8, 5));
}

MATCHER_P(FloatNearPointwise, tol, "Out of range")
{
    return (std::get<0>(arg) > std::get<1>(arg) - tol && std::get<0>(arg) < std::get<1>(arg) + tol);
}

TEST(NdArray1DTest, AccessWholeArray)
{
    auto nx = 20u;
    NdArrayVector1D<> array1d{nx};
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    std::vector<double> numbers(array1d.size());
    for (auto& v : numbers)
    {
        v = dis(gen);
    }

    for (int i = 0; i < static_cast<int>(array1d.size()); ++i)
    {
        array1d(i) = numbers[static_cast<std::vector<double>::size_type>(i)];
    }

    // TODO find better way to compare...

    for (int i = 0; i < static_cast<int>(array1d.size()); ++i)
    {
        EXPECT_DOUBLE_EQ(numbers[static_cast<std::vector<double>::size_type>(i)], array1d(i));
    }
}




TEST(NdArray2DTest, AccessWholeArray)
{
    auto nx = 20u;
    auto ny = 10u;
    NdArrayVector2D<> array2d{nx, ny};
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    std::vector<double> numbers(array2d.size());
    for (auto& v : numbers)
    {
        v = dis(gen);
    }

    std::vector<double>::size_type idx = 0;
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            array2d(i, j) = numbers[idx++];
        }
    }

    // TODO find better way to compare...
    idx = 0;
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            EXPECT_DOUBLE_EQ(numbers[idx++], array2d(i, j));
        }
    }
}



TEST(NdArray3DTest, AccessWholeArray)
{
    auto nx = 20u;
    auto ny = 10u;
    auto nz = 13u;
    NdArrayVector3D<> array3d{nx, ny, nz};
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1.0, 2.0);
    std::vector<double> numbers(array3d.size());
    for (auto& v : numbers)
    {
        v = dis(gen);
    }

    std::vector<double>::size_type idx = 0;
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                array3d(i, j, k) = numbers[idx++];
            }
        }
    }

    // TODO find better way to compare...
    idx = 0;
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                EXPECT_DOUBLE_EQ(numbers[idx++], array3d(i, j, k));
            }
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
