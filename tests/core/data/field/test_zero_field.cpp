#include "core/data/field/zero_field.hpp"
#include "core/numerics/differential_operators.hpp"

#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

#include "gtest/gtest.h"

#include <type_traits>

using namespace PHARE::core;

using Scalar = MHDQuantity::Scalar;

static_assert(has_physicalQuantity_v<ZeroField<FieldMHD<2>>>);
static_assert(std::is_trivially_copyable_v<ZeroField<FieldMHD<2>>>);
static_assert(std::is_trivially_copyable_v<ZeroVecField<VecFieldMHD<2>>>);


TEST(ZeroField, readsZeroForAnyIndex)
{
    ZeroField<FieldMHD<2>> const zero{Scalar::Bx};

    EXPECT_EQ(zero(MeshIndex<2>{3, 4}), 0.);
    EXPECT_EQ(zero(3u, 4u), 0.);
    EXPECT_EQ(zero(Point<std::size_t, 2>{1, 2}), 0.);
}


TEST(ZeroField, projectsToZero)
{
    MeshIndex<2> const index{5, 5};

    EXPECT_EQ((GridLayoutMHD2D::project<Scalar::Bx, Scalar::rho>(ZeroField<FieldMHD<2>>{Scalar::Bx},
                                                                 index)),
              0.);
    EXPECT_EQ((GridLayoutMHD2D::project<Scalar::By, Scalar::Bx>(ZeroField<FieldMHD<2>>{Scalar::By},
                                                                index)),
              0.);
    EXPECT_EQ((GridLayoutMHD2D::project<Scalar::Bz, Scalar::Ex>(ZeroField<FieldMHD<2>>{Scalar::Bz},
                                                                index)),
              0.);
}


TEST(ZeroVecField, componentsCarryTheirQuantities)
{
    ZeroVecField<VecFieldMHD<2>> const zero{MHDQuantity::Vector::B};

    EXPECT_EQ(zero.physicalQuantity(), MHDQuantity::Vector::B);
    EXPECT_EQ(zero(Component::X).physicalQuantity(), Scalar::Bx);
    EXPECT_EQ(zero(Component::Y).physicalQuantity(), Scalar::By);
    EXPECT_EQ(zero(Component::Z).physicalQuantity(), Scalar::Bz);

    EXPECT_EQ(GridLayoutMHD2D::centering(zero(Component::X)),
              GridLayoutMHD2D::centering(Scalar::Bx));
    EXPECT_EQ(GridLayoutMHD2D::centering(zero(Component::Y)),
              GridLayoutMHD2D::centering(Scalar::By));
    EXPECT_EQ(GridLayoutMHD2D::centering(zero(Component::Z)),
              GridLayoutMHD2D::centering(Scalar::Bz));
}


TEST(ZeroVecField, divergenceIsZero)
{
    GridLayoutMHD2D const layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};
    ZeroVecField<VecFieldMHD<2>> const zero{MHDQuantity::Vector::B};

    auto const divergence = div(layout, zero);

    EXPECT_EQ(divergence(MeshIndex<2>{5, 5}), 0.);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
