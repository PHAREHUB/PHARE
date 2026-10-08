#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_boundary_condition_antisymmetric.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_dirichlet.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


TEST_F(FieldBC1D, AntiSymmetricScalarEquivalentToDirichletZero)
{
    Grid1D refGrid{"rho_ref", qty, layout.allocSize(qty)};
    Field1D& refField{*(&refGrid)};
    for (std::uint32_t i = 0; i < refGrid.shape()[0]; ++i)
        refField(i) = (i >= physStart && i <= physEnd) ? interiorValue : ghostSentinel;
    FieldBoundaryConditionDirichlet<Field1D, GridLayout1D> dirichlet{
        0.0, DirichletExtrapolation::Linear};
    dirichlet.apply(refField, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);
    dirichlet.apply(refField, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    FieldBoundaryConditionAntiSymmetric<Field1D, GridLayout1D> antisym;
    antisym.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);
    antisym.apply(field, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    for (std::uint32_t i = 0; i < grid.shape()[0]; ++i)
        EXPECT_DOUBLE_EQ(field(i), refField(i)) << "at index " << i;
}


TEST_F(VecFieldBC1D, AntiSymmetricNormalComponentBxSetToNeumann)
{
    FieldBoundaryConditionAntiSymmetric<VecField1D, GridLayout1D> bc;
    bc.apply(B, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    auto& Bx                  = B[0];
    auto bxQty                = HybridQuantity::Scalar::Bx;
    std::uint32_t bxPhysStart = layout.physicalStartIndex(bxQty, Direction::X);
    std::uint32_t bxPhysEnd   = layout.physicalEndIndex(bxQty, Direction::X);
    EXPECT_DOUBLE_EQ(Bx(bxPhysStart - 1), interiorValue);
    EXPECT_DOUBLE_EQ(Bx(bxPhysEnd + 1), interiorValue);
}

TEST_F(VecFieldBC1D, AntiSymmetricTangentialComponentsByBzSetToDirichletZero)
{
    FieldBoundaryConditionAntiSymmetric<VecField1D, GridLayout1D> bc;
    bc.apply(B, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    auto byQty              = HybridQuantity::Scalar::By;
    std::uint32_t byPhysEnd = layout.physicalEndIndex(byQty, Direction::X);
    for (std::size_t comp : {1u, 2u})
    {
        auto& f = B[comp];
        EXPECT_DOUBLE_EQ(f(0), -interiorValue) << "component " << comp << " lower ghost";
        EXPECT_DOUBLE_EQ(f(byPhysEnd + 1), -interiorValue)
            << "component " << comp << " upper ghost";
    }
}


TEST_F(VecFieldBC2D, AntiSymmetricAtXBoundaries)
{
    FieldBoundaryConditionAntiSymmetric<VecField2D, GridLayout2D> bc;
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, 0.0);

    {
        auto& Bx          = B[0];
        auto bxQty        = HybridQuantity::Scalar::Bx;
        std::uint32_t psx = layout.physicalStartIndex(bxQty, Direction::X);
        std::uint32_t pex = layout.physicalEndIndex(bxQty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(bxQty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(bxQty, Direction::Y);
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
        {
            EXPECT_DOUBLE_EQ(Bx(psx - 1, iy), interiorValue) << "Bx lower ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(Bx(pex + 1, iy), interiorValue) << "Bx upper ghost iy=" << iy;
        }
    }

    for (std::size_t comp : {1u, 2u})
    {
        auto& f           = B[comp];
        auto qty          = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psx = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t pex = layout.physicalEndIndex(qty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(qty, Direction::Y);
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
        {
            EXPECT_DOUBLE_EQ(f(psx - 1, iy), -interiorValue)
                << "comp=" << comp << " lower ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(f(pex + 1, iy), -interiorValue)
                << "comp=" << comp << " upper ghost iy=" << iy;
        }
    }
}

TEST_F(VecFieldBC2D, AntiSymmetricAtYBoundaries)
{
    FieldBoundaryConditionAntiSymmetric<VecField2D, GridLayout2D> bc;
    bc.apply(B, BoundaryLocation::YLower, yLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::YUpper, yUpperGhostCellBox2D(), layout, 0.0);

    {
        auto& By          = B[1];
        auto byQty        = HybridQuantity::Scalar::By;
        std::uint32_t psy = layout.physicalStartIndex(byQty, Direction::Y);
        std::uint32_t pey = layout.physicalEndIndex(byQty, Direction::Y);
        std::uint32_t sx  = layout.physicalStartIndex(byQty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(byQty, Direction::X);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
        {
            EXPECT_DOUBLE_EQ(By(ix, psy - 1), interiorValue) << "By lower ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(By(ix, pey + 1), interiorValue) << "By upper ghost ix=" << ix;
        }
    }

    for (std::size_t comp : {0u, 2u})
    {
        auto& f           = B[comp];
        auto qty          = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psy = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t pey = layout.physicalEndIndex(qty, Direction::Y);
        std::uint32_t sx  = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(qty, Direction::X);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
        {
            EXPECT_DOUBLE_EQ(f(ix, psy - 1), -interiorValue)
                << "comp=" << comp << " lower ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(f(ix, pey + 1), -interiorValue)
                << "comp=" << comp << " upper ghost ix=" << ix;
        }
    }
}


TEST_F(VecFieldBC3D, AntiSymmetricAtZBoundaries)
{
    FieldBoundaryConditionAntiSymmetric<VecField3D, GridLayout3D> bc;
    bc.apply(B, BoundaryLocation::ZLower, zLowerGhostCellBox3D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::ZUpper, zUpperGhostCellBox3D(), layout, 0.0);

    {
        auto& Bz          = B[2];
        auto bzQty        = HybridQuantity::Scalar::Bz;
        std::uint32_t psz = layout.physicalStartIndex(bzQty, Direction::Z);
        std::uint32_t pez = layout.physicalEndIndex(bzQty, Direction::Z);
        std::uint32_t sx  = layout.physicalStartIndex(bzQty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(bzQty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(bzQty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(bzQty, Direction::Y);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iy = sy; iy <= ey; ++iy)
            {
                EXPECT_DOUBLE_EQ(Bz(ix, iy, psz - 1), interiorValue)
                    << "Bz lower ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(Bz(ix, iy, pez + 1), interiorValue)
                    << "Bz upper ghost ix=" << ix << " iy=" << iy;
            }
    }

    for (std::size_t comp : {0u, 1u})
    {
        auto& f           = B[comp];
        auto qty          = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psz = layout.physicalStartIndex(qty, Direction::Z);
        std::uint32_t pez = layout.physicalEndIndex(qty, Direction::Z);
        std::uint32_t sx  = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(qty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(qty, Direction::Y);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iy = sy; iy <= ey; ++iy)
            {
                EXPECT_DOUBLE_EQ(f(ix, iy, psz - 1), -interiorValue)
                    << "comp=" << comp << " lower ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(f(ix, iy, pez + 1), -interiorValue)
                    << "comp=" << comp << " upper ghost ix=" << ix << " iy=" << iy;
            }
    }
}


TEST_F(VecFieldBC3D, AntiSymmetricAtXBoundaries)
{
    FieldBoundaryConditionAntiSymmetric<VecField3D, GridLayout3D> bc;
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox3D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox3D(), layout, 0.0);

    {
        auto& Bx          = B[0];
        auto bxQty        = HybridQuantity::Scalar::Bx;
        std::uint32_t psx = layout.physicalStartIndex(bxQty, Direction::X);
        std::uint32_t pex = layout.physicalEndIndex(bxQty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(bxQty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(bxQty, Direction::Y);
        std::uint32_t sz  = layout.physicalStartIndex(bxQty, Direction::Z);
        std::uint32_t ez  = layout.physicalEndIndex(bxQty, Direction::Z);
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
            for (std::uint32_t iz = sz; iz <= ez; ++iz)
            {
                EXPECT_DOUBLE_EQ(Bx(psx - 1, iy, iz), interiorValue);
                EXPECT_DOUBLE_EQ(Bx(pex + 1, iy, iz), interiorValue);
            }
    }

    for (std::size_t comp : {1u, 2u})
    {
        auto& f           = B[comp];
        auto qty          = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psx = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t pex = layout.physicalEndIndex(qty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(qty, Direction::Y);
        std::uint32_t sz  = layout.physicalStartIndex(qty, Direction::Z);
        std::uint32_t ez  = layout.physicalEndIndex(qty, Direction::Z);
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
            for (std::uint32_t iz = sz; iz <= ez; ++iz)
            {
                EXPECT_DOUBLE_EQ(f(psx - 1, iy, iz), -interiorValue) << "comp=" << comp;
                EXPECT_DOUBLE_EQ(f(pex + 1, iy, iz), -interiorValue) << "comp=" << comp;
            }
    }
}

TEST_F(VecFieldBC3D, AntiSymmetricAtYBoundaries)
{
    FieldBoundaryConditionAntiSymmetric<VecField3D, GridLayout3D> bc;
    bc.apply(B, BoundaryLocation::YLower, yLowerGhostCellBox3D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::YUpper, yUpperGhostCellBox3D(), layout, 0.0);

    {
        auto& By          = B[1];
        auto byQty        = HybridQuantity::Scalar::By;
        std::uint32_t psy = layout.physicalStartIndex(byQty, Direction::Y);
        std::uint32_t pey = layout.physicalEndIndex(byQty, Direction::Y);
        std::uint32_t sx  = layout.physicalStartIndex(byQty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(byQty, Direction::X);
        std::uint32_t sz  = layout.physicalStartIndex(byQty, Direction::Z);
        std::uint32_t ez  = layout.physicalEndIndex(byQty, Direction::Z);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iz = sz; iz <= ez; ++iz)
            {
                EXPECT_DOUBLE_EQ(By(ix, psy - 1, iz), interiorValue);
                EXPECT_DOUBLE_EQ(By(ix, pey + 1, iz), interiorValue);
            }
    }

    for (std::size_t comp : {0u, 2u})
    {
        auto& f           = B[comp];
        auto qty          = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psy = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t pey = layout.physicalEndIndex(qty, Direction::Y);
        std::uint32_t sx  = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(qty, Direction::X);
        std::uint32_t sz  = layout.physicalStartIndex(qty, Direction::Z);
        std::uint32_t ez  = layout.physicalEndIndex(qty, Direction::Z);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iz = sz; iz <= ez; ++iz)
            {
                EXPECT_DOUBLE_EQ(f(ix, psy - 1, iz), -interiorValue) << "comp=" << comp;
                EXPECT_DOUBLE_EQ(f(ix, pey + 1, iz), -interiorValue) << "comp=" << comp;
            }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
