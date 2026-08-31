#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


TEST_F(FieldBC1D, DirichletSetsLowerGhostByLinearExtrapolation)
{
    double const value = 3.0;
    FieldDirichletBoundaryCondition<Field1D, GridLayout1D> bc{value};
    bc.apply(field, BoundaryLocation::XLower, lowerGhostCellBox(), layout, 0.0);

    double expected = 2.0 * value - interiorValue;
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(g), expected);
}

TEST_F(FieldBC1D, DirichletSetsUpperGhostByLinearExtrapolation)
{
    double const value = 3.0;
    FieldDirichletBoundaryCondition<Field1D, GridLayout1D> bc{value};
    bc.apply(field, BoundaryLocation::XUpper, upperGhostCellBox(), layout, 0.0);

    double expected       = 2.0 * value - interiorValue;
    std::uint32_t allocSz = grid.shape()[0];
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        EXPECT_DOUBLE_EQ(field(allocSz - 1 - g), expected);
}


TEST_F(FieldBC2D, DirichletAtXBoundaries)
{
    double const value    = 3.0;
    double const expected = 2.0 * value - interiorValue;
    FieldDirichletBoundaryCondition<Field2D, GridLayout2D> bc{value};
    bc.apply(field, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(field, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, 0.0);

    std::uint32_t const allocX = grid.shape()[0];
    std::uint32_t sy           = layout.physicalStartIndex(qty, Direction::Y);
    std::uint32_t ey           = layout.physicalEndIndex(qty, Direction::Y);
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        for (std::uint32_t iy = sy; iy <= ey; ++iy)
        {
            EXPECT_DOUBLE_EQ(field(g, iy), expected) << "lower g=" << g << " iy=" << iy;
            EXPECT_DOUBLE_EQ(field(allocX - 1 - g, iy), expected) << "upper g=" << g;
        }
}

TEST_F(FieldBC2D, DirichletAtYBoundaries)
{
    double const value    = 3.0;
    double const expected = 2.0 * value - interiorValue;
    FieldDirichletBoundaryCondition<Field2D, GridLayout2D> bc{value};
    bc.apply(field, BoundaryLocation::YLower, yLowerGhostCellBox2D(), layout, 0.0);
    bc.apply(field, BoundaryLocation::YUpper, yUpperGhostCellBox2D(), layout, 0.0);

    std::uint32_t const allocY = grid.shape()[1];
    std::uint32_t sx           = layout.physicalStartIndex(qty, Direction::X);
    std::uint32_t ex           = layout.physicalEndIndex(qty, Direction::X);
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
        {
            EXPECT_DOUBLE_EQ(field(ix, g), expected) << "lower g=" << g << " ix=" << ix;
            EXPECT_DOUBLE_EQ(field(ix, allocY - 1 - g), expected) << "upper g=" << g;
        }
}


TEST_F(FieldBC3D, DirichletAtZBoundaries)
{
    double const value    = 3.0;
    double const expected = 2.0 * value - interiorValue;
    FieldDirichletBoundaryCondition<Field3D, GridLayout3D> bc{value};
    bc.apply(field, BoundaryLocation::ZLower, zLowerGhostCellBox3D(), layout, 0.0);
    bc.apply(field, BoundaryLocation::ZUpper, zUpperGhostCellBox3D(), layout, 0.0);

    std::uint32_t const allocZ = grid.shape()[2];
    std::uint32_t sx           = layout.physicalStartIndex(qty, Direction::X);
    std::uint32_t ex           = layout.physicalEndIndex(qty, Direction::X);
    std::uint32_t sy           = layout.physicalStartIndex(qty, Direction::Y);
    std::uint32_t ey           = layout.physicalEndIndex(qty, Direction::Y);
    for (std::uint32_t g = 0; g < ghostWidth; ++g)
        for (std::uint32_t ix = sx; ix <= ex; ++ix)
            for (std::uint32_t iy = sy; iy <= ey; ++iy)
            {
                EXPECT_DOUBLE_EQ(field(ix, iy, g), expected) << "lower g=" << g;
                EXPECT_DOUBLE_EQ(field(ix, iy, allocZ - 1 - g), expected) << "upper g=" << g;
            }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
