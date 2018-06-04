#include "data/field/field.h"
#include "data/ndarray/ndarray_vector.h"
#include "hybrid/hybrid_quantities.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <ctype.h>
#include <string>

using PHARE::Field;
using PHARE::NdArrayVector1D;
using PHARE::NdArrayVector2D;
using PHARE::NdArrayVector3D;



template<class NdArrayImpl>
class GenericField1D : public ::testing::Test
{
public:
    GenericField1D()
        : f{"test", PHARE::HybridQuantity::Scalar::rho, nx}
    {
    }

protected:
    const uint32_t nx = 10;
    Field<NdArrayImpl, PHARE::HybridQuantity::Scalar> f;
};


template<class NdArrayImpl>
class GenericField2D : public ::testing::Test
{
public:
    GenericField2D()
        : f{"test", PHARE::HybridQuantity::Scalar::rho, nx, ny}
    {
    }

protected:
    const uint32_t nx = 10u;
    const uint32_t ny = 12u;
    Field<NdArrayImpl, PHARE::HybridQuantity::Scalar> f;
};


template<class NdArrayImpl>
class GenericField3D : public ::testing::Test
{
public:
    GenericField3D()
        : f{"test", PHARE::HybridQuantity::Scalar::rho, nx, ny, nz}
    {
    }

protected:
    const uint32_t nx = 10;
    const uint32_t ny = 12;
    const uint32_t nz = 12;
    Field<NdArrayImpl, PHARE::HybridQuantity::Scalar> f;
};


using NdArrays1D = ::testing::Types<NdArrayVector1D<>>;
using NdArrays2D = ::testing::Types<NdArrayVector2D<>>;
using NdArrays3D = ::testing::Types<NdArrayVector3D<>>;

TYPED_TEST_CASE(GenericField1D, NdArrays1D);
TYPED_TEST_CASE(GenericField2D, NdArrays2D);
TYPED_TEST_CASE(GenericField3D, NdArrays3D);



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
    EXPECT_EQ(PHARE::HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}

TYPED_TEST(GenericField2D, physiscalQuantity)
{
    EXPECT_EQ(PHARE::HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}

TYPED_TEST(GenericField3D, physiscalQuantity)
{
    EXPECT_EQ(PHARE::HybridQuantity::Scalar::rho, this->f.physicalQuantity());
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
