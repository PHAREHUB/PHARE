#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

#include <core/data/field/field.h>
#include <core/data/ndarray/ndarray_vector.h>
#include <core/data/vecfield/vecfield.h>
#include <core/hybrid/hybrid_quantities.h>


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
        : vf2{vf2_name, {{HybridQuantity::Bx, HybridQuantity::By, HybridQuantity::Bz}}}
    {
    }

protected:
    std::string vf2_name = "vf";
    VecField<NdArrayImpl, HybridQuantity> vf2;
};

using NdArrays = ::testing::Types<NdArrayVector1D<>, NdArrayVector2D<>, NdArrayVector3D<>>;


TYPED_TEST_CASE(VecFieldGeneric, NdArrays);




class VecFieldTest : public ::testing::Test
{
public:
    VecFieldTest()
        : bx1d_{"field", HybridQuantity::Bx, nx}
        , by1d_{"field", HybridQuantity::By, nx}
        , bz1d_{"field", HybridQuantity::Bz, nx}
        , bx2d_{"field", HybridQuantity::Bx, nx, ny}
        , by2d_{"field", HybridQuantity::By, nx, ny}
        , bz2d_{"field", HybridQuantity::Bz, nx, ny}
        , bx3d_{"field", HybridQuantity::Bz, nx, ny, nz}
        , by3d_{"field", HybridQuantity::Bz, nx, ny, nz}
        , bz3d_{"field", HybridQuantity::Bz, nx, ny, nz}
        , B1D_{"B1D", {{HybridQuantity::Bx, HybridQuantity::By, HybridQuantity::Bz}}}
        , B2D_{"B2D", {{HybridQuantity::Bx, HybridQuantity::By, HybridQuantity::Bz}}}
        , B3D_{"B3D", {{HybridQuantity::Bx, HybridQuantity::By, HybridQuantity::Bz}}}
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

    static const int nx;
    static const int ny;
    static const int nz;
    Field<NdArrayVector1D<>, HybridQuantity> bx1d_;
    Field<NdArrayVector1D<>, HybridQuantity> by1d_;
    Field<NdArrayVector1D<>, HybridQuantity> bz1d_;

    Field<NdArrayVector2D<>, HybridQuantity> bx2d_;
    Field<NdArrayVector2D<>, HybridQuantity> by2d_;
    Field<NdArrayVector2D<>, HybridQuantity> bz2d_;

    Field<NdArrayVector3D<>, HybridQuantity> bx3d_;
    Field<NdArrayVector3D<>, HybridQuantity> by3d_;
    Field<NdArrayVector3D<>, HybridQuantity> bz3d_;

    VecField<NdArrayVector1D<>, HybridQuantity> B1D_;
    VecField<NdArrayVector2D<>, HybridQuantity> B2D_;
    VecField<NdArrayVector3D<>, HybridQuantity> B3D_;
};

const int VecFieldTest::nx = 10;
const int VecFieldTest::ny = 20;
const int VecFieldTest::nz = 30;


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


TEST_F(VecFieldTest, SizeIsOkAfterSet1D)
{
    this->setBuffers();
    EXPECT_EQ(nx, B1D_.getComponent(decltype(B1D_)::Component::X).size());
    EXPECT_EQ(nx, B1D_.getComponent(decltype(B1D_)::Component::Y).size());
    EXPECT_EQ(nx, B1D_.getComponent(decltype(B1D_)::Component::Z).size());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, SizeIsOkAfterSet2D)
{
    this->setBuffers();
    EXPECT_EQ(nx * ny, B2D_.getComponent(decltype(B2D_)::Component::X).size());
    EXPECT_EQ(nx * ny, B2D_.getComponent(decltype(B2D_)::Component::Y).size());
    EXPECT_EQ(nx * ny, B2D_.getComponent(decltype(B2D_)::Component::Z).size());
    this->unsetBuffers();
}


TEST_F(VecFieldTest, SizeIsOkAfterSet3D)
{
    this->setBuffers();
    EXPECT_EQ(nx * ny * nz, B3D_.getComponent(decltype(B3D_)::Component::X).size());
    EXPECT_EQ(nx * ny * nz, B3D_.getComponent(decltype(B3D_)::Component::Y).size());
    EXPECT_EQ(nx * ny * nz, B3D_.getComponent(decltype(B3D_)::Component::Z).size());
    this->unsetBuffers();
}

TEST_F(VecFieldTest, HasThreeBuffers1D)
{
    auto pairs = B1D_.buffersField();
    EXPECT_EQ(3, pairs.size());
}

TEST_F(VecFieldTest, HasThreeBuffers2D)
{
    auto pairs = B2D_.buffersField();
    EXPECT_EQ(3, pairs.size());
}

TEST_F(VecFieldTest, HasThreeBuffers3D)
{
    auto pairs = B3D_.buffersField();
    EXPECT_EQ(3, pairs.size());
}


TEST_F(VecFieldTest, ComponentNames)
{
    auto pairs1D = B1D_.buffersField();
    auto pairs2D = B2D_.buffersField();
    auto pairs3D = B3D_.buffersField();

    EXPECT_EQ(B1D_.name() + "_x", pairs1D[0].first);
    EXPECT_EQ(B1D_.name() + "_y", pairs1D[1].first);
    EXPECT_EQ(B1D_.name() + "_z", pairs1D[2].first);

    EXPECT_EQ(B2D_.name() + "_x", pairs2D[0].first);
    EXPECT_EQ(B2D_.name() + "_y", pairs2D[1].first);
    EXPECT_EQ(B2D_.name() + "_z", pairs2D[2].first);

    EXPECT_EQ(B3D_.name() + "_x", pairs3D[0].first);
    EXPECT_EQ(B3D_.name() + "_y", pairs3D[1].first);
    EXPECT_EQ(B3D_.name() + "_z", pairs3D[2].first);
}



TEST_F(VecFieldTest, PhysicalQuantities)
{
    auto pairs1D = B1D_.buffersField();
    auto pairs2D = B2D_.buffersField();
    auto pairs3D = B3D_.buffersField();

    EXPECT_EQ(HybridQuantity::Bx, pairs1D[0].second);
    EXPECT_EQ(HybridQuantity::By, pairs1D[1].second);
    EXPECT_EQ(HybridQuantity::Bz, pairs1D[2].second);

    EXPECT_EQ(HybridQuantity::Bx, pairs2D[0].second);
    EXPECT_EQ(HybridQuantity::By, pairs2D[1].second);
    EXPECT_EQ(HybridQuantity::Bz, pairs2D[2].second);

    EXPECT_EQ(HybridQuantity::Bx, pairs3D[0].second);
    EXPECT_EQ(HybridQuantity::By, pairs3D[1].second);
    EXPECT_EQ(HybridQuantity::Bz, pairs3D[2].second);
}




int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
