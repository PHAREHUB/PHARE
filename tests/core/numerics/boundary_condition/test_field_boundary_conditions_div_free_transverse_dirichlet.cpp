#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_boundary_condition_divergence_free_transverse_dirichlet.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


TEST_F(VecFieldBC1D, DivergenceFreeTransverseDirichletAtXBoundaries)
{
    std::array values{123.0, 7.0, 11.0};
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecField1D, GridLayout1D> bc{
        values, DirichletExtrapolation::Linear};
    bc.apply(B, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    auto& Bx                  = B[0];
    auto bxQty                = HybridQuantity::Scalar::Bx;
    std::uint32_t bxPhysStart = layout.physicalStartIndex(bxQty, Direction::X);
    std::uint32_t bxPhysEnd   = layout.physicalEndIndex(bxQty, Direction::X);
    EXPECT_DOUBLE_EQ(Bx(bxPhysStart - 2), interiorValue);
    EXPECT_DOUBLE_EQ(Bx(bxPhysStart - 1), interiorValue);
    EXPECT_DOUBLE_EQ(Bx(bxPhysEnd + 1), interiorValue);
    EXPECT_DOUBLE_EQ(Bx(bxPhysEnd + 2), interiorValue);

    for (std::size_t comp : {1u, 2u})
    {
        auto& f              = B[comp];
        auto qty             = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psx    = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t pex    = layout.physicalEndIndex(qty, Direction::X);
        double expectedGhost = 2.0 * values[comp] - interiorValue;
        EXPECT_DOUBLE_EQ(f(psx - 2), expectedGhost) << "component " << comp << " lower far ghost";
        EXPECT_DOUBLE_EQ(f(psx - 1), expectedGhost) << "component " << comp << " lower near ghost";
        EXPECT_DOUBLE_EQ(f(pex + 1), expectedGhost) << "component " << comp << " upper near ghost";
        EXPECT_DOUBLE_EQ(f(pex + 2), expectedGhost) << "component " << comp << " upper far ghost";
    }
}


TEST_F(VecFieldBC2D, DivergenceFreeTransverseDirichletAtXBoundaries)
{
    std::array values{123.0, 7.0, 11.0};
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecField2D, GridLayout2D> bc{
        values, DirichletExtrapolation::Linear};
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, 0.0);

    auto& Bx          = B[0];
    auto bxQty        = HybridQuantity::Scalar::Bx;
    std::uint32_t psx = layout.physicalStartIndex(bxQty, Direction::X);
    std::uint32_t pex = layout.physicalEndIndex(bxQty, Direction::X);
    std::uint32_t sy  = layout.physicalStartIndex(bxQty, Direction::Y);
    std::uint32_t ey  = layout.physicalEndIndex(bxQty, Direction::Y);
    for (std::uint32_t iy = sy; iy <= ey; ++iy)
    {
        EXPECT_DOUBLE_EQ(Bx(psx - 2, iy), interiorValue) << "Bx lower far ghost iy=" << iy;
        EXPECT_DOUBLE_EQ(Bx(psx - 1, iy), interiorValue) << "Bx lower near ghost iy=" << iy;
        EXPECT_DOUBLE_EQ(Bx(pex + 1, iy), interiorValue) << "Bx upper near ghost iy=" << iy;
        EXPECT_DOUBLE_EQ(Bx(pex + 2, iy), interiorValue) << "Bx upper far ghost iy=" << iy;
    }

    for (std::size_t comp : {1u, 2u})
    {
        auto& f              = B[comp];
        auto qty             = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psx    = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t pex    = layout.physicalEndIndex(qty, Direction::X);
        std::uint32_t sy     = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t ey     = layout.physicalEndIndex(qty, Direction::Y);
        double expectedGhost = 2.0 * values[comp] - interiorValue;
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
        {
            EXPECT_DOUBLE_EQ(f(psx - 2, iy), expectedGhost)
                << "comp=" << comp << " lower far ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(f(psx - 1, iy), expectedGhost)
                << "comp=" << comp << " lower near ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(f(pex + 1, iy), expectedGhost)
                << "comp=" << comp << " upper near ghost iy=" << iy;
            EXPECT_DOUBLE_EQ(f(pex + 2, iy), expectedGhost)
                << "comp=" << comp << " upper far ghost iy=" << iy;
        }
    }
}

TEST_F(VecFieldBC2D, DivergenceFreeTransverseDirichletAtYBoundaries)
{
    std::array values{3.0, 123.0, 11.0};
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecField2D, GridLayout2D> bc{
        values, DirichletExtrapolation::Linear};
    bc.apply(B, BoundaryLocation::YLower, yLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::YUpper, yUpperGhostCellBox2D(), layout, 0.0);

    auto& By          = B[1];
    auto byQty        = HybridQuantity::Scalar::By;
    std::uint32_t psy = layout.physicalStartIndex(byQty, Direction::Y);
    std::uint32_t pey = layout.physicalEndIndex(byQty, Direction::Y);
    std::uint32_t sx  = layout.physicalStartIndex(byQty, Direction::X);
    std::uint32_t ex  = layout.physicalEndIndex(byQty, Direction::X);
    for (std::uint32_t ix = sx; ix <= ex; ++ix)
    {
        EXPECT_DOUBLE_EQ(By(ix, psy - 2), interiorValue) << "By lower far ghost ix=" << ix;
        EXPECT_DOUBLE_EQ(By(ix, psy - 1), interiorValue) << "By lower near ghost ix=" << ix;
        EXPECT_DOUBLE_EQ(By(ix, pey + 1), interiorValue) << "By upper near ghost ix=" << ix;
        EXPECT_DOUBLE_EQ(By(ix, pey + 2), interiorValue) << "By upper far ghost ix=" << ix;
    }

    for (std::size_t comp : {0u, 2u})
    {
        auto& f              = B[comp];
        auto qty             = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psy    = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t pey    = layout.physicalEndIndex(qty, Direction::Y);
        std::uint32_t sx     = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t ex     = layout.physicalEndIndex(qty, Direction::X);
        double expectedGhost = 2.0 * values[comp] - interiorValue;
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
        {
            EXPECT_DOUBLE_EQ(f(ix, psy - 2), expectedGhost)
                << "comp=" << comp << " lower far ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(f(ix, psy - 1), expectedGhost)
                << "comp=" << comp << " lower near ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(f(ix, pey + 1), expectedGhost)
                << "comp=" << comp << " upper near ghost ix=" << ix;
            EXPECT_DOUBLE_EQ(f(ix, pey + 2), expectedGhost)
                << "comp=" << comp << " upper far ghost ix=" << ix;
        }
    }
}

TEST_F(VecFieldBC2DNonUniformBy, DivergenceFreeTransverseDirichletKeepsXGhostDivergenceZero)
{
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecField2D, GridLayout2D> bc{
        std::array{123.0, 0.0, 11.0}, DirichletExtrapolation::Linear};
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, 0.0);

    auto& Bx = B[0];
    auto& By = B[1];
    for (auto const& index : xLowerGhostCellBox2D())
    {
        EXPECT_NEAR(Bx(index.shift(0, 1)) - Bx(index) + By(index.shift(1, 1)) - By(index), 0.0,
                    1e-12)
            << "lower divergence at (" << index[0] << ", " << index[1] << ")";
    }
    for (auto const& index : xUpperGhostCellBox2D())
    {
        EXPECT_NEAR(Bx(index.shift(0, 1)) - Bx(index) + By(index.shift(1, 1)) - By(index), 0.0,
                    1e-12)
            << "upper divergence at (" << index[0] << ", " << index[1] << ")";
    }
}


TEST_F(VecFieldBC2DNonUniformByAnisotropic,
       DivergenceFreeTransverseDirichletKeepsXGhostDivergenceZeroOnAnisotropicMesh)
{
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecField2D, GridLayout2D> bc{
        std::array{123.0, 0.0, 11.0}, DirichletExtrapolation::Linear};
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, 0.0);

    auto& Bx        = B[0];
    auto& By        = B[1];
    double const dx = layout.meshSize()[0];
    double const dy = layout.meshSize()[1];
    auto divB       = [&](auto const& index) {
        return (Bx(index.shift(0, 1)) - Bx(index)) / dx + (By(index.shift(1, 1)) - By(index)) / dy;
    };
    for (auto const& index : xLowerGhostCellBox2D())
        EXPECT_NEAR(divB(index), 0.0, 1e-12)
            << "lower divergence at (" << index[0] << ", " << index[1] << ")";
    for (auto const& index : xUpperGhostCellBox2D())
        EXPECT_NEAR(divB(index), 0.0, 1e-12)
            << "upper divergence at (" << index[0] << ", " << index[1] << ")";
}


TEST_F(VecFieldBC3D, DivergenceFreeTransverseDirichletAtZBoundaries)
{
    std::array values{3.0, 7.0, 123.0};
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecField3D, GridLayout3D> bc{
        values, DirichletExtrapolation::Linear};
    bc.apply(B, BoundaryLocation::ZLower, zLowerGhostCellBox3D(), layout, 0.0);
    bc.apply(B, BoundaryLocation::ZUpper, zUpperGhostCellBox3D(), layout, 0.0);

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
            EXPECT_DOUBLE_EQ(Bz(ix, iy, psz - 2), interiorValue)
                << "Bz lower far ghost ix=" << ix << " iy=" << iy;
            EXPECT_DOUBLE_EQ(Bz(ix, iy, psz - 1), interiorValue)
                << "Bz lower near ghost ix=" << ix << " iy=" << iy;
            EXPECT_DOUBLE_EQ(Bz(ix, iy, pez + 1), interiorValue)
                << "Bz upper near ghost ix=" << ix << " iy=" << iy;
            EXPECT_DOUBLE_EQ(Bz(ix, iy, pez + 2), interiorValue)
                << "Bz upper far ghost ix=" << ix << " iy=" << iy;
        }

    for (std::size_t comp : {0u, 1u})
    {
        auto& f              = B[comp];
        auto qty             = HybridQuantity::componentsQuantities(vecQty)[comp];
        std::uint32_t psz    = layout.physicalStartIndex(qty, Direction::Z);
        std::uint32_t pez    = layout.physicalEndIndex(qty, Direction::Z);
        std::uint32_t sx     = layout.physicalStartIndex(qty, Direction::X);
        std::uint32_t ex     = layout.physicalEndIndex(qty, Direction::X);
        std::uint32_t sy     = layout.physicalStartIndex(qty, Direction::Y);
        std::uint32_t ey     = layout.physicalEndIndex(qty, Direction::Y);
        double expectedGhost = 2.0 * values[comp] - interiorValue;
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iy = sy; iy <= ey; ++iy)
            {
                EXPECT_DOUBLE_EQ(f(ix, iy, psz - 2), expectedGhost)
                    << "comp=" << comp << " lower far ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(f(ix, iy, psz - 1), expectedGhost)
                    << "comp=" << comp << " lower near ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(f(ix, iy, pez + 1), expectedGhost)
                    << "comp=" << comp << " upper near ghost ix=" << ix << " iy=" << iy;
                EXPECT_DOUBLE_EQ(f(ix, iy, pez + 2), expectedGhost)
                    << "comp=" << comp << " upper far ghost ix=" << ix << " iy=" << iy;
            }
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
