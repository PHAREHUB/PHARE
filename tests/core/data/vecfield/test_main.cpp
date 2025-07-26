#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <string>

#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/hybrid/hybrid_quantities.hpp"


using namespace PHARE::core;



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
    VecField<Field<NdArrayImpl::dimension, HybridQuantity::Scalar>, HybridQuantity> vf2;
};

using NdArrays = ::testing::Types<NdArrayView<1>, NdArrayView<2>, NdArrayView<3>>;

TYPED_TEST_SUITE(VecFieldGeneric, NdArrays);

template<std::size_t dim>
using Field_t = Field<dim, HybridQuantity::Scalar>;

template<std::size_t dim>
using VecField_t = VecField<Field_t<dim>, HybridQuantity>;


class VecFieldTest : public ::testing::Test
{
public:
    VecFieldTest()
        : bx1d_{"B1D_x", HybridQuantity::Scalar::Bx, nx}
        , by1d_{"B1D_y", HybridQuantity::Scalar::By, nx}
        , bz1d_{"B1D_z", HybridQuantity::Scalar::Bz, nx}
        , bx2d_{"B2D_x", HybridQuantity::Scalar::Bx, nx, ny}
        , by2d_{"B2D_y", HybridQuantity::Scalar::By, nx, ny}
        , bz2d_{"B2D_z", HybridQuantity::Scalar::Bz, nx, ny}
        , bx3d_{"B3D_x", HybridQuantity::Scalar::Bx, nx, ny, nz}
        , by3d_{"B3D_y", HybridQuantity::Scalar::By, nx, ny, nz}
        , bz3d_{"B3D_z", HybridQuantity::Scalar::Bz, nx, ny, nz}
        , B1D_{"B1D", HybridQuantity::Vector::B}
        , B2D_{"B2D", HybridQuantity::Vector::B}
        , B3D_{"B3D", HybridQuantity::Vector::B}
    {
    }

protected:
    void setBuffers()
    {
        B1D_[0].setBuffer(&bx1d_);
        B1D_[1].setBuffer(&by1d_);
        B1D_[2].setBuffer(&bz1d_);

        B2D_[0].setBuffer(&bx2d_);
        B2D_[1].setBuffer(&by2d_);
        B2D_[2].setBuffer(&bz2d_);

        B3D_[0].setBuffer(&bx3d_);
        B3D_[1].setBuffer(&by3d_);
        B3D_[2].setBuffer(&bz3d_);
    }

    void unsetBuffers()
    {
        B1D_[0].setBuffer(nullptr);
        B1D_[1].setBuffer(nullptr);
        B1D_[2].setBuffer(nullptr);

        B2D_[0].setBuffer(nullptr);
        B2D_[1].setBuffer(nullptr);
        B2D_[2].setBuffer(nullptr);

        B3D_[0].setBuffer(nullptr);
        B3D_[1].setBuffer(nullptr);
        B3D_[2].setBuffer(nullptr);
    }

    static const std::uint32_t nx;
    static const std::uint32_t ny;
    static const std::uint32_t nz;
    Grid<NdArrayVector<1>, typename HybridQuantity::Scalar> bx1d_;
    Grid<NdArrayVector<1>, typename HybridQuantity::Scalar> by1d_;
    Grid<NdArrayVector<1>, typename HybridQuantity::Scalar> bz1d_;

    Grid<NdArrayVector<2>, typename HybridQuantity::Scalar> bx2d_;
    Grid<NdArrayVector<2>, typename HybridQuantity::Scalar> by2d_;
    Grid<NdArrayVector<2>, typename HybridQuantity::Scalar> bz2d_;

    Grid<NdArrayVector<3>, typename HybridQuantity::Scalar> bx3d_;
    Grid<NdArrayVector<3>, typename HybridQuantity::Scalar> by3d_;
    Grid<NdArrayVector<3>, typename HybridQuantity::Scalar> bz3d_;

    VecField_t<1> B1D_;
    VecField_t<2> B2D_;
    VecField_t<3> B3D_;
};

const std::uint32_t VecFieldTest::nx = 10;
const std::uint32_t VecFieldTest::ny = 20;
const std::uint32_t VecFieldTest::nz = 30;


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




TEST_F(VecFieldTest, VecFieldsHaveBeginAndEnd)
{
    for ([[maybe_unused]] auto const& component : B1D_) {}
}




TEST(aVecField, dataCanBeCopiedIntoAnother)
{
    using Scalar = typename HybridQuantity::Scalar;

    Grid<NdArrayVector<3>, Scalar> bx1{"B1_x", Scalar::Bx, 2u, 3u, 4u};
    Grid<NdArrayVector<3>, Scalar> by1{"B1_y", Scalar::By, 2u, 3u, 4u};
    Grid<NdArrayVector<3>, Scalar> bz1{"B1_z", Scalar::Bz, 2u, 3u, 4u};
    VecField_t<3> B1{"B1", HybridQuantity::Vector::B};
    B1[0].setBuffer(&bx1);
    B1[1].setBuffer(&by1);
    B1[2].setBuffer(&bz1);

    bx1(1, 1, 1) = 12;
    by1(1, 1, 1) = 13;
    bz1(1, 1, 1) = 14;

    Grid<NdArrayVector<3>, Scalar> bx2{"B2_x", Scalar::Bx, 2u, 3u, 4u};
    Grid<NdArrayVector<3>, Scalar> by2{"B2_y", Scalar::By, 2u, 3u, 4u};
    Grid<NdArrayVector<3>, Scalar> bz2{"B2_z", Scalar::Bz, 2u, 3u, 4u};
    VecField_t<3> B2{"B2", HybridQuantity::Vector::B};
    B2[0].setBuffer(&bx2);
    B2[1].setBuffer(&by2);
    B2[2].setBuffer(&bz2);

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
