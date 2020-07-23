
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <random>
#include <string>

#include "core/data/ndarray/ndarray_vector.h"


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
    const std::uint32_t nx = 10;
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
    const std::uint32_t nx = 10;
    const std::uint32_t ny = 20;
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
    const std::uint32_t nx = 10;
    const std::uint32_t ny = 20;
    const std::uint32_t nz = 30;
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



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
