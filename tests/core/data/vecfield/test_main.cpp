#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

#include "data/field/field.h"
#include "data/ndarray/ndarray_vector.h"
#include "data/vecfield/vecfield.h"
#include "hybrid/hybrid_quantities.h"


using PHARE::Component;
using PHARE::Field;
using PHARE::HybridQuantity;
using PHARE::NdArrayVector1D;
using PHARE::NdArrayVector2D;
using PHARE::NdArrayVector3D;
using PHARE::VecField;



template<class NdArrayImpl>
class VecFieldGeneric : public ::testing::Test
{
public:
    VecFieldGeneric()
        : vf2{vf2_name,
              {{HybridQuantity::Scalar::Bx, HybridQuantity::Scalar::By,
                HybridQuantity::Scalar::Bz}}}
    {
    }

protected:
    std::string vf2_name = "vf";
    VecField<NdArrayImpl, HybridQuantity::Scalar> vf2;
};

using NdArrays = ::testing::Types<NdArrayVector1D<>, NdArrayVector2D<>, NdArrayVector3D<>>;


TYPED_TEST_CASE(VecFieldGeneric, NdArrays);




class VecFieldTest : public ::testing::Test
{
public:
    VecFieldTest()
        : bx1d_{"field", HybridQuantity::Scalar::Bx, nx}
        , by1d_{"field", HybridQuantity::Scalar::By, nx}
        , bz1d_{"field", HybridQuantity::Scalar::Bz, nx}
        , bx2d_{"field", HybridQuantity::Scalar::Bx, nx, ny}
        , by2d_{"field", HybridQuantity::Scalar::By, nx, ny}
        , bz2d_{"field", HybridQuantity::Scalar::Bz, nx, ny}
        , bx3d_{"field", HybridQuantity::Scalar::Bx, nx, ny, nz}
        , by3d_{"field", HybridQuantity::Scalar::By, nx, ny, nz}
        , bz3d_{"field", HybridQuantity::Scalar::Bz, nx, ny, nz}
        , B1D_{"B1D", HybridQuantity::Vector::B}
        , B2D_{"B2D", HybridQuantity::Vector::B}
        , B3D_{"B3D", HybridQuantity::Vector::B}
    {
    }

protected:
    void setBuffers()
    {
        B1D_.setBuffer("B1D_x", &bx1d_);
        B1D_.setBuffer("B1D_y", &by1d_);
        B1D_.setBuffer("B1D_z", &bz1d_);

        B2D_.setBuffer("B2D_x", &bx2d_);
        B2D_.setBuffer("B2D_y", &by2d_);
        B2D_.setBuffer("B2D_z", &bz2d_);

        B3D_.setBuffer("B3D_x", &bx3d_);
        B3D_.setBuffer("B3D_y", &by3d_);
        B3D_.setBuffer("B3D_z", &bz3d_);
    }

    void unsetBuffers()
    {
        B1D_.setBuffer("B1D_x", nullptr);
        B1D_.setBuffer("B1D_y", nullptr);
        B1D_.setBuffer("B1D_z", nullptr);

        B2D_.setBuffer("B2D_x", nullptr);
        B2D_.setBuffer("B2D_y", nullptr);
        B2D_.setBuffer("B2D_z", nullptr);

        B3D_.setBuffer("B3D_x", nullptr);
        B3D_.setBuffer("B3D_y", nullptr);
        B3D_.setBuffer("B3D_z", nullptr);
    }

    static const uint32_t nx;
    static const uint32_t ny;
    static const uint32_t nz;
    Field<NdArrayVector1D<>, typename HybridQuantity::Scalar> bx1d_;
    Field<NdArrayVector1D<>, typename HybridQuantity::Scalar> by1d_;
    Field<NdArrayVector1D<>, typename HybridQuantity::Scalar> bz1d_;

    Field<NdArrayVector2D<>, typename HybridQuantity::Scalar> bx2d_;
    Field<NdArrayVector2D<>, typename HybridQuantity::Scalar> by2d_;
    Field<NdArrayVector2D<>, typename HybridQuantity::Scalar> bz2d_;

    Field<NdArrayVector3D<>, typename HybridQuantity::Scalar> bx3d_;
    Field<NdArrayVector3D<>, typename HybridQuantity::Scalar> by3d_;
    Field<NdArrayVector3D<>, typename HybridQuantity::Scalar> bz3d_;

    VecField<NdArrayVector1D<>, HybridQuantity> B1D_;
    VecField<NdArrayVector2D<>, HybridQuantity> B2D_;
    VecField<NdArrayVector3D<>, HybridQuantity> B3D_;
};

const uint32_t VecFieldTest::nx = 10;
const uint32_t VecFieldTest::ny = 20;
const uint32_t VecFieldTest::nz = 30;


TEST_F(VecFieldTest, isNotInitiallyUsable1D)
{
    EXPECT_FALSE(B1D_.isUsable());
}


TEST_F(VecFieldTest, isNotInitiallyUsable2D)
{
    EXPECT_FALSE(B2D_.isUsable());
}


TEST_F(VecFieldTest, isNotInitiallyUsable3D)
{
    EXPECT_FALSE(B3D_.isUsable());
}


TEST_F(VecFieldTest, isSettableInitially1D)
{
    EXPECT_TRUE(B1D_.isSettable());
}

TEST_F(VecFieldTest, isSettableInitially2D)
{
    EXPECT_TRUE(B2D_.isSettable());
}


TEST_F(VecFieldTest, isSettableInitially3D)
{
    EXPECT_TRUE(B3D_.isSettable());
}



TEST_F(VecFieldTest, IsUsableAfterSet1D)
{
    this->setBuffers();
    EXPECT_TRUE(B1D_.isUsable());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, IsUsableAfterSet2D)
{
    this->setBuffers();
    EXPECT_TRUE(B2D_.isUsable());
    this->unsetBuffers();
}



TEST_F(VecFieldTest, IsUsableAfterSet3D)
{
    this->setBuffers();
    EXPECT_TRUE(B3D_.isUsable());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, isNotSettableAfterSet1D)
{
    this->setBuffers();
    EXPECT_FALSE(B1D_.isSettable());
    this->unsetBuffers();
}

TEST_F(VecFieldTest, isNotSettableAfterSet2D)
{
    this->setBuffers();
    EXPECT_FALSE(B2D_.isSettable());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, isNotSettableAfterSet3D)
{
    this->setBuffers();
    EXPECT_FALSE(B3D_.isSettable());
    this->unsetBuffers();
}




TEST_F(VecFieldTest, SizeIsOkAfterSet1D)
{
    this->setBuffers();
    EXPECT_EQ(nx, B1D_.getComponent(Component::X).size());
    EXPECT_EQ(nx, B1D_.getComponent(Component::Y).size());
    EXPECT_EQ(nx, B1D_.getComponent(Component::Z).size());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, SizeIsOkAfterSet2D)
{
    this->setBuffers();
    EXPECT_EQ(nx * ny, B2D_.getComponent(Component::X).size());
    EXPECT_EQ(nx * ny, B2D_.getComponent(Component::Y).size());
    EXPECT_EQ(nx * ny, B2D_.getComponent(Component::Z).size());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, SizeIsOkAfterSet3D)
{
    this->setBuffers();
    EXPECT_EQ(nx * ny * nz, B3D_.getComponent(Component::X).size());
    EXPECT_EQ(nx * ny * nz, B3D_.getComponent(Component::Y).size());
    EXPECT_EQ(nx * ny * nz, B3D_.getComponent(Component::Z).size());
    this->unsetBuffers();
}

TEST_F(VecFieldTest, HasThreeBuffers1D)
{
    auto pairs = B1D_.getFieldNamesAndQuantities();
    EXPECT_EQ(3, pairs.size());
}

TEST_F(VecFieldTest, HasThreeBuffers2D)
{
    auto pairs = B2D_.getFieldNamesAndQuantities();
    EXPECT_EQ(3, pairs.size());
}

TEST_F(VecFieldTest, HasThreeBuffers3D)
{
    auto pairs = B3D_.getFieldNamesAndQuantities();
    EXPECT_EQ(3, pairs.size());
}


TEST_F(VecFieldTest, ComponentNames)
{
    auto pairs1D = B1D_.getFieldNamesAndQuantities();
    auto pairs2D = B2D_.getFieldNamesAndQuantities();
    auto pairs3D = B3D_.getFieldNamesAndQuantities();

    EXPECT_EQ(B1D_.name() + "_x", pairs1D[0].name);
    EXPECT_EQ(B1D_.name() + "_y", pairs1D[1].name);
    EXPECT_EQ(B1D_.name() + "_z", pairs1D[2].name);

    EXPECT_EQ(B2D_.name() + "_x", pairs2D[0].name);
    EXPECT_EQ(B2D_.name() + "_y", pairs2D[1].name);
    EXPECT_EQ(B2D_.name() + "_z", pairs2D[2].name);

    EXPECT_EQ(B3D_.name() + "_x", pairs3D[0].name);
    EXPECT_EQ(B3D_.name() + "_y", pairs3D[1].name);
    EXPECT_EQ(B3D_.name() + "_z", pairs3D[2].name);
}



TEST_F(VecFieldTest, PhysicalQuantities)
{
    auto pairs1D = B1D_.getFieldNamesAndQuantities();
    auto pairs2D = B2D_.getFieldNamesAndQuantities();
    auto pairs3D = B3D_.getFieldNamesAndQuantities();

    EXPECT_EQ(HybridQuantity::Scalar::Bx, pairs1D[0].qty);
    EXPECT_EQ(HybridQuantity::Scalar::By, pairs1D[1].qty);
    EXPECT_EQ(HybridQuantity::Scalar::Bz, pairs1D[2].qty);

    EXPECT_EQ(HybridQuantity::Scalar::Bx, pairs2D[0].qty);
    EXPECT_EQ(HybridQuantity::Scalar::By, pairs2D[1].qty);
    EXPECT_EQ(HybridQuantity::Scalar::Bz, pairs2D[2].qty);

    EXPECT_EQ(HybridQuantity::Scalar::Bx, pairs3D[0].qty);
    EXPECT_EQ(HybridQuantity::Scalar::By, pairs3D[1].qty);
    EXPECT_EQ(HybridQuantity::Scalar::Bz, pairs3D[2].qty);
}



TEST(aVecField, dataCanBeCopiedIntoAnother)
{
    using Scalar = typename HybridQuantity::Scalar;

    Field<NdArrayVector3D<>, Scalar> bx1{"B1_bx1", Scalar::Bx, 2u, 3u, 4u};
    Field<NdArrayVector3D<>, Scalar> by1{"B1_by1", Scalar::By, 2u, 3u, 4u};
    Field<NdArrayVector3D<>, Scalar> bz1{"B1_bz1", Scalar::Bz, 2u, 3u, 4u};
    VecField<NdArrayVector3D<>, HybridQuantity> B1{"B1", HybridQuantity::Vector::B};
    B1.setBuffer("B1_x", &bx1);
    B1.setBuffer("B1_y", &by1);
    B1.setBuffer("B1_z", &bz1);

    bx1(1, 1, 1) = 12;
    by1(1, 1, 1) = 13;
    bz1(1, 1, 1) = 14;

    Field<NdArrayVector3D<>, Scalar> bx2{"B2_bx2", Scalar::Bx, 2u, 3u, 4u};
    Field<NdArrayVector3D<>, Scalar> by2{"B2_by2", Scalar::By, 2u, 3u, 4u};
    Field<NdArrayVector3D<>, Scalar> bz2{"B2_bz2", Scalar::Bz, 2u, 3u, 4u};
    VecField<NdArrayVector3D<>, HybridQuantity> B2{"B2", HybridQuantity::Vector::B};
    B2.setBuffer("B2_x", &bx2);
    B2.setBuffer("B2_y", &by2);
    B2.setBuffer("B2_z", &bz2);

    B2.copyData(B1);

    EXPECT_DOUBLE_EQ(12, bx2(1, 1, 1));
    EXPECT_DOUBLE_EQ(13, by2(1, 1, 1));
    EXPECT_DOUBLE_EQ(14, bz2(1, 1, 1));
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
