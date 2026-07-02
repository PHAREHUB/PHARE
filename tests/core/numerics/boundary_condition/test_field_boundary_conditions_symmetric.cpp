#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_symmetric_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


// ─── 1D scalar ────────────────────────────────────────────────────────────────

TEST_F(FieldBC1D, SymmetricScalarEquivalentToNeumann)
{
    // Reference: apply Neumann on a copy
    Grid1D refGrid{"rho_ref", qty, layout.allocSize(qty)};
    Field1D& refField{*(&refGrid)};
    for (std::uint32_t i = 0; i < refGrid.shape()[0]; ++i)
        refField(i) = (i >= physStart && i <= physEnd) ? interiorValue : ghostSentinel;
    FieldNeumannBoundaryCondition<Field1D, GridLayout1D> neumann;
    neumann.apply(refField, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    neumann.apply(refField, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

    FieldSymmetricBoundaryCondition<Field1D, GridLayout1D> sym;
    sym.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    sym.apply(field, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

    for (std::uint32_t i = 0; i < grid.shape()[0]; ++i)
        EXPECT_DOUBLE_EQ(field(i), refField(i)) << "at index " << i;
}


// ─── 1D VecField ─────────────────────────────────────────────────────────────

TEST_F(VecFieldBC1D, SymmetricNormalComponentBxSetToDirichletZero)
{
    FieldSymmetricBoundaryCondition<VecField1D, GridLayout1D> bc;
    bc.apply(B, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto& Bx                  = B[0];
    auto bxQty                = HybridQuantity::Scalar::Bx;
    std::uint32_t bxPhysStart = layout.physicalStartIndex(bxQty, Direction::X);
    std::uint32_t bxPhysEnd   = layout.physicalEndIndex(bxQty, Direction::X);
    EXPECT_DOUBLE_EQ(Bx(bxPhysStart - 1), -interiorValue);
    EXPECT_DOUBLE_EQ(Bx(bxPhysStart), 0.0);
    EXPECT_DOUBLE_EQ(Bx(bxPhysEnd), 0.0);
    EXPECT_DOUBLE_EQ(Bx(bxPhysEnd + 1), -interiorValue);
}

TEST_F(VecFieldBC1D, SymmetricTangentialComponentsByBzSetToNeumann)
{
    FieldSymmetricBoundaryCondition<VecField1D, GridLayout1D> bc;
    bc.apply(B, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto byQty              = HybridQuantity::Scalar::By;
    std::uint32_t byPhysEnd = layout.physicalEndIndex(byQty, Direction::X);
    for (std::size_t comp : {1u, 2u})
    {
        auto& f = B[comp];
        EXPECT_DOUBLE_EQ(f(0), interiorValue) << "component " << comp << " lower ghost";
        EXPECT_DOUBLE_EQ(f(byPhysEnd + 1), interiorValue) << "component " << comp << " upper ghost";
    }
}


// ─── 2D VecField ─────────────────────────────────────────────────────────────

TEST_F(VecFieldBC2D, SymmetricAtXBoundaries)
{
    FieldSymmetricBoundaryCondition<VecField2D, GridLayout2D> bc;
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, makeCtx(acc, 0.0));

    // Bx: primal in X → Dirichlet(0)
    {
        auto& Bx          = B[0];
        auto bxQty        = HybridQuantity::Scalar::Bx;
        std::uint32_t psx = layout.physicalStartIndex(bxQty, Direction::X);
        std::uint32_t pex = layout.physicalEndIndex(bxQty, Direction::X);
        std::uint32_t sy  = layout.physicalStartIndex(bxQty, Direction::Y);
        std::uint32_t ey  = layout.physicalEndIndex(bxQty, Direction::Y);
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
        {
            EXPECT_DOUBLE_EQ(Bx(psx - 1, iy), -interiorValue) << "Bx lower ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(Bx(psx, iy), 0.0) << "Bx lower boundary iy=" << iy;
            EXPECT_DOUBLE_EQ(Bx(pex, iy), 0.0) << "Bx upper boundary iy=" << iy;
            EXPECT_DOUBLE_EQ(Bx(pex + 1, iy), -interiorValue) << "Bx upper ghost iy=" << iy;
        }
    }

    // By, Bz: dual in X → Neumann
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
            EXPECT_DOUBLE_EQ(f(psx - 1, iy), interiorValue)
                << "comp=" << comp << " lower ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(f(pex + 1, iy), interiorValue)
                << "comp=" << comp << " upper ghost iy=" << iy;
        }
    }
}

TEST_F(VecFieldBC2D, SymmetricAtYBoundaries)
{
    FieldSymmetricBoundaryCondition<VecField2D, GridLayout2D> bc;
    bc.apply(B, BoundaryLocation::YLower, yLowerGhostCellBox2D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::YUpper, yUpperGhostCellBox2D(), layout, makeCtx(acc, 0.0));

    // By: primal in Y → Dirichlet(0)
    {
        auto& By          = B[1];
        auto byQty        = HybridQuantity::Scalar::By;
        std::uint32_t psy = layout.physicalStartIndex(byQty, Direction::Y);
        std::uint32_t pey = layout.physicalEndIndex(byQty, Direction::Y);
        std::uint32_t sx  = layout.physicalStartIndex(byQty, Direction::X);
        std::uint32_t ex  = layout.physicalEndIndex(byQty, Direction::X);
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
        {
            EXPECT_DOUBLE_EQ(By(ix, psy - 1), -interiorValue) << "By lower ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(By(ix, psy), 0.0) << "By lower boundary ix=" << ix;
            EXPECT_DOUBLE_EQ(By(ix, pey), 0.0) << "By upper boundary ix=" << ix;
            EXPECT_DOUBLE_EQ(By(ix, pey + 1), -interiorValue) << "By upper ghost ix=" << ix;
        }
    }

    // Bx, Bz: dual in Y → Neumann
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
            EXPECT_DOUBLE_EQ(f(ix, psy - 1), interiorValue)
                << "comp=" << comp << " lower ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(f(ix, pey + 1), interiorValue)
                << "comp=" << comp << " upper ghost ix=" << ix;
        }
    }
}


// ─── 3D VecField ─────────────────────────────────────────────────────────────

TEST_F(VecFieldBC3D, SymmetricAtZBoundaries)
{
    FieldSymmetricBoundaryCondition<VecField3D, GridLayout3D> bc;
    bc.apply(B, BoundaryLocation::ZLower, zLowerGhostCellBox3D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::ZUpper, zUpperGhostCellBox3D(), layout, makeCtx(acc, 0.0));

    // Bz: primal in Z → Dirichlet(0)
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
                EXPECT_DOUBLE_EQ(Bz(ix, iy, psz - 1), -interiorValue)
                    << "Bz lower ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(Bz(ix, iy, psz), 0.0)
                    << "Bz lower boundary ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(Bz(ix, iy, pez), 0.0)
                    << "Bz upper boundary ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(Bz(ix, iy, pez + 1), -interiorValue)
                    << "Bz upper ghost ix=" << ix << " iy=" << iy;
            }
    }

    // Bx, By: dual in Z → Neumann
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
                EXPECT_DOUBLE_EQ(f(ix, iy, psz - 1), interiorValue)
                    << "comp=" << comp << " lower ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(f(ix, iy, pez + 1), interiorValue)
                    << "comp=" << comp << " upper ghost ix=" << ix << " iy=" << iy;
            }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
