#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <random>
#include <string>

#include "core/data/ndarray/ndarray_vector.hpp"


using namespace PHARE::core;


template<class NdArray>
class GenericNdArray1D : public ::testing::Test
{
public:
    GenericNdArray1D()
        : a{{nx}}
    {
    }

protected:
    std::uint32_t const nx = 10;
    NdArray a;
};


template<class NdArray>
class GenericNdArray2D : public ::testing::Test
{
public:
    GenericNdArray2D()
        : a{{nx, ny}}
    {
    }

protected:
    std::uint32_t const nx = 10;
    std::uint32_t const ny = 20;
    NdArray a;
};


template<class NdArray>
class GenericNdArray3D : public ::testing::Test
{
public:
    GenericNdArray3D()
        : a{{nx, ny, nz}}
    {
    }

protected:
    std::uint32_t const nx = 10;
    std::uint32_t const ny = 20;
    std::uint32_t const nz = 30;
    NdArray a;
};


using NdArray1D = ::testing::Types<NdArrayVector<1>>;
using NdArray2D = ::testing::Types<NdArrayVector<2>>;
using NdArray3D = ::testing::Types<NdArrayVector<3>>;


TYPED_TEST_SUITE(GenericNdArray1D, NdArray1D);
TYPED_TEST_SUITE(GenericNdArray2D, NdArray2D);
TYPED_TEST_SUITE(GenericNdArray3D, NdArray3D);



TYPED_TEST(GenericNdArray1D, SizeIsOkAfterNonEmptyCreation)
{
    EXPECT_EQ(this->nx, this->a.size());
}


TYPED_TEST(GenericNdArray1D, IsModifiable)
{
    std::uint32_t i{2};
    this->a(i) = 12.;
    EXPECT_EQ(12., this->a(i));
}


TYPED_TEST(GenericNdArray1D, CanBeReadOnly)
{
    std::uint32_t i{2};
    this->a(i)                  = 12.;
    NdArrayVector<1> const& ref = this->a;
    EXPECT_EQ(12., ref(i));
}


TYPED_TEST(GenericNdArray1D, AccessWholeArray)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0, 2.0);
    std::vector<double> numbers(this->a.size());

    for (auto& v : numbers)
        v = dis(gen);

    for (std::uint32_t i = 0; i < this->nx; ++i)
    {
        this->a(i) = numbers[i];
        EXPECT_DOUBLE_EQ(numbers[i], this->a(i));
    }
}


TYPED_TEST(GenericNdArray1D, HasCopyCtor)
{
    for (auto& e : this->a)
        e = 12.;

    NdArrayVector<1> other{this->nx};
    other = this->a;

    for (auto const& e : other)
        EXPECT_DOUBLE_EQ(12., e);
}


TYPED_TEST(GenericNdArray1D, HasMoveCtor)
{
    for (auto& e : this->a)
        e = 12.;

    NdArrayVector<1> other = std::move(this->a);

    for (auto const& e : other)
        EXPECT_DOUBLE_EQ(12., e);
}



TYPED_TEST(GenericNdArray2D, SizeIsOkAfterNonEmptyCreation)
{
    EXPECT_EQ(this->nx * this->ny, this->a.size());
}


TYPED_TEST(GenericNdArray2D, IsModifiable)
{
    std::uint32_t i{2}, j{3};
    this->a(i, j) = 12.;
    EXPECT_EQ(12., this->a(i, j));
}


TYPED_TEST(GenericNdArray2D, CanBeReadOnly)
{
    std::uint32_t i{2}, j{3};
    this->a(i, j)               = 12.;
    NdArrayVector<2> const& ref = this->a;
    EXPECT_EQ(12., ref(i, j));
}


TYPED_TEST(GenericNdArray2D, AccessWholeArray)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0, 2.0);
    std::vector<double> numbers(this->a.size());

    for (auto& v : numbers)
        v = dis(gen);

    std::uint32_t idx{0};
    for (std::uint32_t i = 0; i < this->nx; ++i)
    {
        for (std::uint32_t j = 0; j < this->ny; ++j)
        {
            this->a(i, j) = numbers[idx];
            EXPECT_DOUBLE_EQ(numbers[idx++], this->a(i, j));
        }
    }
}


TYPED_TEST(GenericNdArray2D, HasCopyCtor)
{
    for (auto& e : this->a)
        e = 12.;

    NdArrayVector<2> other{this->nx, this->ny};
    other = this->a;

    for (auto const& e : other)
        EXPECT_DOUBLE_EQ(12., e);
}


TYPED_TEST(GenericNdArray2D, HasMoveCtor)
{
    for (auto& e : this->a)
        e = 12.;

    NdArrayVector<2> other = std::move(this->a);

    for (auto const& e : other)
        EXPECT_DOUBLE_EQ(12., e);
}



TYPED_TEST(GenericNdArray3D, SizeIsOkAfterNonEmptyCreation)
{
    EXPECT_EQ(this->nx * this->ny * this->nz, this->a.size());
}


TYPED_TEST(GenericNdArray3D, IsModifiable)
{
    std::uint32_t i{2}, j{3}, k{4};
    this->a(i, j, k) = 12.;
    EXPECT_EQ(12., this->a(i, j, k));
}


TYPED_TEST(GenericNdArray3D, CanBeReadOnly)
{
    std::uint32_t i{2}, j{3}, k{4};
    this->a(i, j, k)            = 12.;
    NdArrayVector<3> const& ref = this->a;
    EXPECT_EQ(12., ref(i, j, k));
}


TYPED_TEST(GenericNdArray3D, AccessWholeArray)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0, 2.0);
    std::vector<double> numbers(this->a.size());

    for (auto& v : numbers)
        v = dis(gen);


    std::uint32_t idx{0};
    for (std::uint32_t i = 0; i < this->nx; ++i)
    {
        for (std::uint32_t j = 0; j < this->ny; ++j)
        {
            for (std::uint32_t k = 0; k < this->nz; ++k)
            {
                this->a(i, j, k) = numbers[idx];
                EXPECT_DOUBLE_EQ(numbers[idx++], this->a(i, j, k));
            }
        }
    }
}


TYPED_TEST(GenericNdArray3D, HasCopyCtor)
{
    for (auto& e : this->a)
        e = 12.;

    NdArrayVector<3> other{this->nx, this->ny, this->nz};
    other = this->a;

    for (auto const& e : other)
        EXPECT_DOUBLE_EQ(12., e);
}


TYPED_TEST(GenericNdArray3D, HasMoveCtor)
{
    for (auto& e : this->a)
        e = 12.;

    NdArrayVector<3> other = std::move(this->a);

    for (auto const& e : other)
        EXPECT_DOUBLE_EQ(12., e);
}



TEST(MaskedView1d, maskOps)
{
    constexpr std::size_t dim    = 1;
    constexpr std::uint32_t size = 20;
    using Mask                   = NdArrayMask;
    NdArrayVector<dim> array{{size}, 0.};

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), 0);

    array[Mask{0u, size - 1}] = 1;

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), size);

    Mask oneCellOffset2{2u};
    array[oneCellOffset2] = 2;

    EXPECT_EQ(2, oneCellOffset2.nCells(array));
    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), size + oneCellOffset2.nCells(array));

    array[Mask{5u, 6u}] = 2;

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), size + 6);

    EXPECT_EQ(array(0), 1);
    EXPECT_EQ(array(size - 1), 1);
    array[Mask{5u}] >> array[Mask{0u}];
    EXPECT_EQ(array(0), 2);
    EXPECT_EQ(array(size - 1), 2);

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), size + 8);
}

TEST(MaskedView2d, maskOps)
{
    constexpr std::size_t dim      = 2;
    constexpr std::uint32_t size   = 20;
    constexpr std::uint32_t sizeSq = 20 * 20;
    using Mask                     = NdArrayMask;
    NdArrayVector<dim> array{{size, size}, 0.};

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), 0);

    std::fill(array.begin(), array.end(), 1);

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), sizeSq);

    Mask oneCellOffset2{2u};
    array[oneCellOffset2] = 2;

    EXPECT_EQ(60, oneCellOffset2.nCells(array));
    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0),
              sizeSq + oneCellOffset2.nCells(array));

    Mask twoCellsOffset5{5u, 6u};
    array[twoCellsOffset5] = 2;

    EXPECT_EQ((8 * 4 + 4) + (6 * 4 + 4), twoCellsOffset5.nCells(array));
    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0),
              sizeSq + oneCellOffset2.nCells(array) + twoCellsOffset5.nCells(array));

    EXPECT_EQ(array(0, 0), 1);
    EXPECT_EQ(array(size - 1, size - 1), 1);
    array[Mask{5u}] >> array[Mask{0u}];
    EXPECT_EQ(array(0, 0), 2);
    EXPECT_EQ(array(size - 1, size - 1), 2);

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), sizeSq + oneCellOffset2.nCells(array)
                                                                  + twoCellsOffset5.nCells(array)
                                                                  + Mask{0u}.nCells(array));
}

TEST(MaskedView2d, maskOps2)
{
    constexpr std::size_t dim     = 2;
    constexpr std::uint32_t size0 = 20, size1 = 22;
    constexpr std::uint32_t sizeSq = size0 * size1;
    using Mask                     = NdArrayMask;
    NdArrayVector<dim> array{{size0, size1}, 0.};

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), 0);

    std::fill(array.begin(), array.end(), 1);

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), sizeSq);

    Mask oneCellOffset2{2u};
    array[oneCellOffset2] = 2;

    EXPECT_EQ(14 * 2 + 16 * 2 + 4, oneCellOffset2.nCells(array));
    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0),
              sizeSq + oneCellOffset2.nCells(array));

    Mask twoCellsOffset5{5u, 6u};
    array[twoCellsOffset5] = 2;

    EXPECT_EQ((8 * 2 + 10 * 2 + 4) + (6 * 2 + 8 * 2 + 4), twoCellsOffset5.nCells(array));
    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0),
              sizeSq + oneCellOffset2.nCells(array) + twoCellsOffset5.nCells(array));

    EXPECT_EQ(array(0, 0), 1);
    EXPECT_EQ(array(size0 - 1, size1 - 1), 1);
    array[Mask{5u}] >> array[Mask{0u}];
    EXPECT_EQ(array(0, 0), 2);
    EXPECT_EQ(array(size0 - 1, size1 - 1), 2);

    EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), sizeSq + oneCellOffset2.nCells(array)
                                                                  + twoCellsOffset5.nCells(array)
                                                                  + Mask{0u}.nCells(array));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
