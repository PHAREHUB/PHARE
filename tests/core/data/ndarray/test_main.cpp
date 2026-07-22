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

TEST(MaskedView3d, maskOps3)
{
    constexpr std::size_t dim      = 3;
    constexpr std::uint32_t size0  = 10;
    constexpr std::uint32_t sizeCu = size0 * size0 * size0;
    using Mask                     = PHARE::core::NdArrayMask;


    {
        NdArrayVector<dim> array{{size0, size0, size0}, 0.};
        EXPECT_EQ(std::accumulate(array.begin(), array.end(), 0), 0);
        array[Mask{0}] = 1;

        // outter cells of a 10**3 cube = 10**3 - 8**3 = 488
        EXPECT_EQ(sum(array), 488);


        std::fill(array.begin(), array.end(), 1);
        array[Mask{0}] = 0;
        EXPECT_EQ(sum(array), sizeCu - 488);
        array[Mask{1}] >> array[Mask{0}];
        EXPECT_EQ(sum(array), sizeCu);
    }


    PHARE::core::NdArrayVector<3> array({size0, size0, size0}, 0.);

    array[Mask{0}] = 1;
    EXPECT_EQ(sum(array), 488);
    array[Mask{1}] >> array[Mask{0}];
    EXPECT_EQ(sum(array), 0);

    array[Mask{2}] = 1;
    EXPECT_EQ(sum(array), 152);
    array[Mask{1}] = 1;
    EXPECT_EQ(sum(array), 448);
    array[Mask{1}] = 0;
    EXPECT_EQ(sum(array), 152);

    array[Mask{2}] >> array[Mask{1}];
    EXPECT_EQ(sum(array), 448);
    array[Mask{2}] = 0;
    EXPECT_EQ(sum(array), 296);

    EXPECT_EQ(Mask{1}.nCells(array), 296);
    EXPECT_EQ(Mask{2}.nCells(array), 152);
}

TEST(MaskedView3d, operatorRightShiftCoordinatesCorrect)
{
    constexpr std::size_t dim     = 3;
    constexpr std::uint32_t size0 = 10;
    using Mask                    = PHARE::core::NdArrayMask;

    // Create source and destination arrays
    NdArrayVector<dim> outer{{size0, size0, size0}, 0.};
    NdArrayVector<dim> inner{{size0, size0, size0}, 0.};

    // Fill outer's Mask{1} region (second ghost layer) with 1.0
    outer[Mask{1}] = 1.0;

    // Copy from outer's Mask{1} to inner's Mask{0}
    // This exercises all 6 face sections of operator>>
    outer[Mask{1}] >> inner[Mask{0}];

    // If coordinates are correct, all 488 cells of Mask{0} should be 1.0

    EXPECT_DOUBLE_EQ(sum(inner), 488.0);

    // Detailed cell-by-cell verification to pinpoint failures
    for (std::uint32_t i = 0; i < size0; ++i)
    {
        for (std::uint32_t j = 0; j < size0; ++j)
        {
            for (std::uint32_t k = 0; k < size0; ++k)
            {
                bool is_mask0 = (i == 0 || i == size0 - 1 || j == 0 || j == size0 - 1 || k == 0
                                 || k == size0 - 1);
                if (is_mask0)
                {
                    EXPECT_DOUBLE_EQ(inner(i, j, k), 1.0)
                        << "Mask{0} cell (" << i << "," << j << "," << k
                        << ") should be 1.0 but is " << inner(i, j, k);
                }
            }
        }
    }
}


TEST(NdArrayViewFillFrom, sameOrderingCopiesAllValues)
{
    constexpr std::uint32_t nx = 2, ny = 3;
    std::vector<double> src(nx * ny), dst(nx * ny, 0.);
    NdArrayView<2, double> srcv{src.data(), {nx, ny}};
    NdArrayView<2, double> dstv{dst.data(), {nx, ny}};

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            srcv(i, j) = 10. * i + j;

    dstv.fill_from(srcv);

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            EXPECT_DOUBLE_EQ(dstv(i, j), srcv(i, j));
}


TEST(NdArrayViewFillFrom, oppositeOrdering2dTransposesLayout)
{
    // non-square shape with distinct values: a flat buffer copy would scramble
    // the elements, so this only passes if fill_from goes through both views'
    // own index arithmetic.
    constexpr std::uint32_t nx = 2, ny = 3;
    std::vector<double> fbuf(nx * ny), cbuf(nx * ny, 0.);
    NdArrayView<2, double, /*c_ordering=*/false> fortran{fbuf.data(), {nx, ny}};
    NdArrayView<2, double> c{cbuf.data(), {nx, ny}};

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            fortran(i, j) = 10. * i + j;

    c.fill_from(fortran);

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            EXPECT_DOUBLE_EQ(c(i, j), fortran(i, j));

    // the two flat buffers must actually differ, or this test could not
    // distinguish an index-aware copy from a flat one
    EXPECT_NE(fbuf, cbuf);
}


TEST(NdArrayViewFillFrom, oppositeOrdering2dReverseDirection)
{
    constexpr std::uint32_t nx = 2, ny = 3;
    std::vector<double> cbuf(nx * ny), fbuf(nx * ny, 0.);
    NdArrayView<2, double> c{cbuf.data(), {nx, ny}};
    NdArrayView<2, double, false> fortran{fbuf.data(), {nx, ny}};

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            c(i, j) = 10. * i + j;

    fortran.fill_from(c);

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            EXPECT_DOUBLE_EQ(fortran(i, j), c(i, j));

    EXPECT_NE(fbuf, cbuf);
}


TEST(NdArrayViewFillFrom, oppositeOrdering3dCopiesAllValues)
{
    constexpr std::uint32_t nx = 2, ny = 3, nz = 4;
    std::vector<double> fbuf(nx * ny * nz), cbuf(nx * ny * nz, 0.);
    NdArrayView<3, double, false> fortran{fbuf.data(), {nx, ny, nz}};
    NdArrayView<3, double> c{cbuf.data(), {nx, ny, nz}};

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            for (std::uint32_t k = 0; k < nz; ++k)
                fortran(i, j, k) = 100. * i + 10. * j + k;

    c.fill_from(fortran);

    for (std::uint32_t i = 0; i < nx; ++i)
        for (std::uint32_t j = 0; j < ny; ++j)
            for (std::uint32_t k = 0; k < nz; ++k)
                EXPECT_DOUBLE_EQ(c(i, j, k), fortran(i, j, k));

    EXPECT_NE(fbuf, cbuf);
}


TEST(NdArrayViewFillFrom, oppositeOrdering1dIsExactCopy)
{
    constexpr std::uint32_t nx = 5;
    std::vector<double> fbuf(nx), cbuf(nx, 0.);
    NdArrayView<1, double, false> fortran{fbuf.data(), {nx}};
    NdArrayView<1, double> c{cbuf.data(), {nx}};

    for (std::uint32_t i = 0; i < nx; ++i)
        fortran(i) = static_cast<double>(i);

    c.fill_from(fortran);

    EXPECT_EQ(fbuf, cbuf); // orderings coincide in 1D
}


TEST(NdArrayViewFillFrom, oppositeOrderingShapeMismatchThrows)
{
    std::vector<double> fbuf(2 * 3), cbuf(3 * 2);
    NdArrayView<2, double, false> fortran{fbuf.data(), {2, 3}};
    NdArrayView<2, double> c{cbuf.data(), {3, 2}};

    EXPECT_THROW(c.fill_from(fortran), std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
