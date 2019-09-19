
#include <ctype.h>
#include <string>

#include "data/field/field.h"
#include "data/ndarray/ndarray_vector.h"
#include "hybrid/hybrid_quantities.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;


template<class NdArrayImpl>
class GenericField1D : public ::testing::Test
{
public:
    GenericField1D()
        : f{"test", HybridQuantity::Scalar::rho, nx}
    {
    }

protected:
    const uint32_t nx = 10;
    Field<NdArrayImpl, HybridQuantity::Scalar> f;
};


template<class NdArrayImpl>
class GenericField2D : public ::testing::Test
{
public:
    GenericField2D()
        : f{"test", HybridQuantity::Scalar::rho, nx, ny}
    {
    }

protected:
    const uint32_t nx = 10u;
    const uint32_t ny = 12u;
    Field<NdArrayImpl, HybridQuantity::Scalar> f;
};


template<class NdArrayImpl>
class GenericField3D : public ::testing::Test
{
public:
    GenericField3D()
        : f{"test", HybridQuantity::Scalar::rho, nx, ny, nz}
    {
    }

protected:
    const uint32_t nx = 10;
    const uint32_t ny = 12;
    const uint32_t nz = 12;
    Field<NdArrayImpl, HybridQuantity::Scalar> f;
};


using NdArrays1D = ::testing::Types<NdArrayVector1D<>>;
using NdArrays2D = ::testing::Types<NdArrayVector2D<>>;
using NdArrays3D = ::testing::Types<NdArrayVector3D<>>;

TYPED_TEST_SUITE(GenericField1D, NdArrays1D);
TYPED_TEST_SUITE(GenericField2D, NdArrays2D);
TYPED_TEST_SUITE(GenericField3D, NdArrays3D);



TYPED_TEST(GenericField1D, ReturnsTotalNumberOfElements)
{
    EXPECT_EQ(this->nx, this->f.size());
}

TYPED_TEST(GenericField2D, ReturnsTotalNumberOfElements)
{
    EXPECT_EQ(this->nx * this->ny, this->f.size());
}


TYPED_TEST(GenericField3D, ReturnsTotalNumberOfElements)
{
    EXPECT_EQ(this->nx * this->ny * this->nz, this->f.size());
}


TYPED_TEST(GenericField1D, IsModifiable)
{
    this->f(3) = 4.8;
    EXPECT_FLOAT_EQ(4.8, this->f(3));
}


TYPED_TEST(GenericField2D, IsModifiable)
{
    this->f(3, 8) = 8.4;
    EXPECT_FLOAT_EQ(8.4, this->f(3, 8));
}


TYPED_TEST(GenericField3D, IsModifiable)
{
    this->f(3, 8, 2) = 84.48;
    EXPECT_FLOAT_EQ(84.48, this->f(3, 8, 2));
}



TYPED_TEST(GenericField1D, CanBeReadOnly)
{
    this->f(3)      = 4.8;
    auto const& ref = this->f;
    EXPECT_FLOAT_EQ(4.8, ref(3));
}


TYPED_TEST(GenericField2D, CanBeReadOnly)
{
    this->f(3, 8)   = 8.4;
    auto const& ref = this->f;
    EXPECT_FLOAT_EQ(8.4, ref(3, 8));
}


TYPED_TEST(GenericField3D, CanBeReadOnly)
{
    this->f(3, 8, 2) = 84.48;
    auto const& ref  = this->f;
    EXPECT_FLOAT_EQ(84.48, ref(3, 8, 2));
}




TYPED_TEST(GenericField1D, returnName)
{
    EXPECT_EQ(std::string("test"), this->f.name());
}

TYPED_TEST(GenericField2D, returnName)
{
    EXPECT_EQ(std::string("test"), this->f.name());
}

TYPED_TEST(GenericField3D, returnName)
{
    EXPECT_EQ(std::string("test"), this->f.name());
}



TYPED_TEST(GenericField1D, physiscalQuantity)
{
    EXPECT_EQ(HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}

TYPED_TEST(GenericField2D, physiscalQuantity)
{
    EXPECT_EQ(HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}

TYPED_TEST(GenericField3D, physiscalQuantity)
{
    EXPECT_EQ(HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}


TEST(Field1D, canBeAssigned)
{
    auto nx = 10u;
    Field<NdArrayVector1D<>, HybridQuantity::Scalar> f{"test", HybridQuantity::Scalar::rho, nx};
    Field<NdArrayVector1D<>, HybridQuantity::Scalar> other{"other", HybridQuantity::Scalar::rho,
                                                           nx};

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




TEST(Field2D, canBeAssigned)
{
    auto nx = 10u;
    auto ny = 11u;
    Field<NdArrayVector2D<>, HybridQuantity::Scalar> f{"test", HybridQuantity::Scalar::rho, nx, ny};
    Field<NdArrayVector2D<>, HybridQuantity::Scalar> other{"other", HybridQuantity::Scalar::rho, nx,
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




TEST(Field3D, canBeAssigned)
{
    auto nx = 10u;
    auto ny = 11u;
    auto nz = 12u;
    Field<NdArrayVector3D<>, HybridQuantity::Scalar> f{"test", HybridQuantity::Scalar::rho, nx, ny,
                                                       nz};
    Field<NdArrayVector3D<>, HybridQuantity::Scalar> other{"other", HybridQuantity::Scalar::rho, nx,
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




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
