

#include "core/data/grid/grid.hpp"
#include <core/utilities/algorithm.hpp>
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"


#include "gtest/gtest.h"


#include <string>


using namespace PHARE::core;


template<class NdArrayImpl>
class GenericGrid1D : public ::testing::Test
{
public:
    GenericGrid1D()
        : f{"test", HybridQuantity::Scalar::rho, nx}
    {
    }

protected:
    std::uint32_t const nx = 10;
    Grid<NdArrayImpl, HybridQuantity::Scalar> f;
};


template<class NdArrayImpl>
class GenericGrid2D : public ::testing::Test
{
public:
    GenericGrid2D()
        : f{"test", HybridQuantity::Scalar::rho, nx, ny}
    {
    }

protected:
    std::uint32_t const nx = 10u;
    std::uint32_t const ny = 12u;
    Grid<NdArrayImpl, HybridQuantity::Scalar> f;
};


template<class NdArrayImpl>
class GenericGrid3D : public ::testing::Test
{
public:
    GenericGrid3D()
        : f{"test", HybridQuantity::Scalar::rho, nx, ny, nz}
    {
    }

protected:
    std::uint32_t const nx = 10;
    std::uint32_t const ny = 12;
    std::uint32_t const nz = 12;
    Grid<NdArrayImpl, HybridQuantity::Scalar> f;
};


using NdArrays1D = ::testing::Types<NdArrayVector<1>>;
using NdArrays2D = ::testing::Types<NdArrayVector<2>>;
using NdArrays3D = ::testing::Types<NdArrayVector<3>>;

TYPED_TEST_SUITE(GenericGrid1D, NdArrays1D);
TYPED_TEST_SUITE(GenericGrid2D, NdArrays2D);
TYPED_TEST_SUITE(GenericGrid3D, NdArrays3D);



TYPED_TEST(GenericGrid1D, ReturnsTotalNumberOfElements)
{
    EXPECT_EQ(this->nx, this->f.size());
}

TYPED_TEST(GenericGrid2D, ReturnsTotalNumberOfElements)
{
    EXPECT_EQ(this->nx * this->ny, this->f.size());
}


TYPED_TEST(GenericGrid3D, ReturnsTotalNumberOfElements)
{
    EXPECT_EQ(this->nx * this->ny * this->nz, this->f.size());
}


TYPED_TEST(GenericGrid1D, IsModifiable)
{
    this->f(3u) = 4.8;
    EXPECT_FLOAT_EQ(4.8, this->f(3u));
}


TYPED_TEST(GenericGrid2D, IsModifiable)
{
    this->f(3u, 8u) = 8.4;
    EXPECT_FLOAT_EQ(8.4, this->f(3u, 8u));
}


TYPED_TEST(GenericGrid3D, IsModifiable)
{
    this->f(3u, 8u, 2u) = 84.48;
    EXPECT_FLOAT_EQ(84.48, this->f(3u, 8u, 2u));
}



TYPED_TEST(GenericGrid1D, CanBeReadOnly)
{
    this->f(3u)     = 4.8;
    auto const& ref = this->f;
    EXPECT_FLOAT_EQ(4.8, ref(3u));
}


TYPED_TEST(GenericGrid2D, CanBeReadOnly)
{
    this->f(3u, 8u) = 8.4;
    auto const& ref = this->f;
    EXPECT_FLOAT_EQ(8.4, ref(3u, 8u));
}


TYPED_TEST(GenericGrid3D, CanBeReadOnly)
{
    this->f(3u, 8u, 2u) = 84.48;
    auto const& ref     = this->f;
    EXPECT_FLOAT_EQ(84.48, ref(3u, 8u, 2u));
}




TYPED_TEST(GenericGrid1D, returnName)
{
    EXPECT_EQ(std::string("test"), this->f.name());
}

TYPED_TEST(GenericGrid2D, returnName)
{
    EXPECT_EQ(std::string("test"), this->f.name());
}

TYPED_TEST(GenericGrid3D, returnName)
{
    EXPECT_EQ(std::string("test"), this->f.name());
}



TYPED_TEST(GenericGrid1D, physicalQuantity)
{
    EXPECT_EQ(HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}

TYPED_TEST(GenericGrid2D, physicalQuantity)
{
    EXPECT_EQ(HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}

TYPED_TEST(GenericGrid3D, physicalQuantity)
{
    EXPECT_EQ(HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}




TEST(Grid1D, canBeAssigned)
{
    auto nx = 10u;
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> f{"test", HybridQuantity::Scalar::rho, nx};
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> other{"other", HybridQuantity::Scalar::rho, nx};

    for (auto& v : f)
    {
        v = 12.;
    }

    other.copyData(f);

    for (auto const& v : f)
    {
        EXPECT_DOUBLE_EQ(12., v);
    }
}




TEST(Grid2D, canBeAssigned)
{
    auto nx = 10u;
    auto ny = 11u;
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> f{"test", HybridQuantity::Scalar::rho, nx, ny};
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> other{"other", HybridQuantity::Scalar::rho, nx,
                                                         ny};

    for (auto& v : f)
    {
        v = 12.;
    }

    other.copyData(f);

    for (auto const& v : f)
    {
        EXPECT_DOUBLE_EQ(12., v);
    }
}




TEST(Grid3D, canBeAssigned)
{
    auto nx = 10u;
    auto ny = 11u;
    auto nz = 12u;
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> f{"test", HybridQuantity::Scalar::rho, nx, ny,
                                                     nz};
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> other{"other", HybridQuantity::Scalar::rho, nx,
                                                         ny, nz};

    for (auto& v : f)
    {
        v = 12.;
    }

    other.copyData(f);

    for (auto const& v : f)
    {
        EXPECT_DOUBLE_EQ(12., v);
    }
}




TEST(Grid1D, canBeAveraged)
{
    auto nx = 15u;
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> f1{"f1", HybridQuantity::Scalar::rho, nx};
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> f2{"f2", HybridQuantity::Scalar::rho, nx};
    Grid<NdArrayVector<1>, HybridQuantity::Scalar> avg{"f2", HybridQuantity::Scalar::rho, nx};

    //
    for (auto& v : f1)
    {
        v = 10.;
    }


    for (auto& v : f2)
    {
        v = 20.;
    }

    PHARE::core::average(f1, f2, avg);

    for (auto const& v : avg)
    {
        EXPECT_DOUBLE_EQ(15., v);
    }
}




TEST(Grid2D, canBeAveraged)
{
    auto nx = 15u;
    auto ny = 25u;
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> f1{"f1", HybridQuantity::Scalar::rho, nx, ny};
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> f2{"f2", HybridQuantity::Scalar::rho, nx, ny};
    Grid<NdArrayVector<2>, HybridQuantity::Scalar> avg{"f2", HybridQuantity::Scalar::rho, nx, ny};

    //
    for (auto& v : f1)
    {
        v = 10.;
    }


    for (auto& v : f2)
    {
        v = 20.;
    }

    PHARE::core::average(f1, f2, avg);

    for (auto const& v : avg)
    {
        EXPECT_DOUBLE_EQ(15., v);
    }
}




TEST(Grid3D, canBeAveraged)
{
    auto nx = 15u;
    auto ny = 25u;
    auto nz = 35u;
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> f1{"f1", HybridQuantity::Scalar::rho, nx, ny,
                                                      nz};
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> f2{"f2", HybridQuantity::Scalar::rho, nx, ny,
                                                      nz};
    Grid<NdArrayVector<3>, HybridQuantity::Scalar> avg{"f2", HybridQuantity::Scalar::rho, nx, ny,
                                                       nz};

    //
    for (auto& v : f1)
    {
        v = 10.;
    }


    for (auto& v : f2)
    {
        v = 20.;
    }

    PHARE::core::average(f1, f2, avg);

    for (auto const& v : avg)
    {
        EXPECT_DOUBLE_EQ(15., v);
    }
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
