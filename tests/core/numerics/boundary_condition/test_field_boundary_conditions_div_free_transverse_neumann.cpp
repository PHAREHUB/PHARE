#include "gtest/gtest.h"

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_neumann_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

using namespace PHARE::core;


struct DivFreeTransverseNeumannBC2D : testing::Test
{
    GridLayoutMHD2D layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};


    UsableVecFieldMHD<2> Bvec{"bc_test_B", layout, MHDQuantity::Vector::B};


    DivFreeTransverseNeumannBC2D()
    {
        auto fill = [&](auto& f, double a, double b, double c) {
            auto const shape = f.shape();
            for (std::uint32_t i = 0; i < shape[0]; ++i)
                for (std::uint32_t j = 0; j < shape[1]; ++j)
                    f(i, j) = a * static_cast<double>(i) + b * static_cast<double>(j) + c;
        };
        fill(Bvec[0], 0.3, 0.7, 1.0);
        fill(Bvec[1], -0.2, 0.4, 2.0);
        fill(Bvec[2], 0.5, -0.6, 3.0);
    }

    void checkTransverseMirrored(BoundaryLocation loc, Box<std::uint32_t, 2> const& ghostBox)
    {
        Direction const direction = getDirection(loc);
        Side const side           = getSide(loc);
        std::size_t const iNormal = static_cast<std::size_t>(direction);

        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            if (comp == iNormal)
                continue;
            auto& bField         = Bvec[comp];
            auto const qty       = MHDQuantity::componentsQuantities(MHDQuantity::Vector::B)[comp];
            auto const centering = layout.centering(qty)[iNormal];
            for (auto const& index : layout.toFieldBox(ghostBox, qty))
            {
                auto const mirror = layout.boundaryMirrored(direction, side, centering, index);
                EXPECT_DOUBLE_EQ(bField(index), bField(mirror))
                    << "comp=" << comp << " index=(" << index[0] << "," << index[1] << ")";
            }
        }
    }

    void checkGhostCellsDivergenceFree(BoundaryLocation loc, Box<std::uint32_t, 2> const& ghostBox)
    {
        auto& Bx = Bvec[0];
        auto& By = Bvec[1];
        for (auto const& cell : ghostBox)
        {
            double const div
                = (Bx(cell.shift(0, 1)) - Bx(cell)) + (By(cell.shift(1, 1)) - By(cell));
            EXPECT_NEAR(div, 0.0, 1e-12)
                << "cell=(" << cell[0] << "," << cell[1] << ") at " << static_cast<int>(loc);
        }
    }

    void applyAndCheck(BoundaryLocation loc, Box<std::uint32_t, 2> const& ghostBox)
    {
        auto B = Bvec.super();
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D> bc;
        bc.apply(B, loc, ghostBox, layout, 0.0);

        checkTransverseMirrored(loc, ghostBox);
        checkGhostCellsDivergenceFree(loc, ghostBox);
    }
};

TEST_F(DivFreeTransverseNeumannBC2D, XLowerGhostsMirrorTransverseAndKillDivergence)
{
    applyAndCheck(BoundaryLocation::XLower, mhd2DXLowerGhostBox());
}

TEST_F(DivFreeTransverseNeumannBC2D, XUpperGhostsMirrorTransverseAndKillDivergence)
{
    applyAndCheck(BoundaryLocation::XUpper, mhd2DXUpperGhostBox());
}

TEST_F(DivFreeTransverseNeumannBC2D, YLowerGhostsMirrorTransverseAndKillDivergence)
{
    applyAndCheck(BoundaryLocation::YLower, mhd2DYLowerGhostBox());
}

TEST_F(DivFreeTransverseNeumannBC2D, YUpperGhostsMirrorTransverseAndKillDivergence)
{
    applyAndCheck(BoundaryLocation::YUpper, mhd2DYUpperGhostBox());
}


struct DivFreeTransverseNeumannBC2DAnisotropic : testing::Test
{
    GridLayoutMHD2D layout{{0.1, 0.2}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};


    UsableVecFieldMHD<2> Bvec{"bc_test_B", layout, MHDQuantity::Vector::B};


    DivFreeTransverseNeumannBC2DAnisotropic()
    {
        auto fill = [&](auto& f, double a, double b, double c) {
            auto const shape = f.shape();
            for (std::uint32_t i = 0; i < shape[0]; ++i)
                for (std::uint32_t j = 0; j < shape[1]; ++j)
                    f(i, j) = a * static_cast<double>(i) + b * static_cast<double>(j) + c;
        };
        fill(Bvec[0], 0.3, 0.7, 1.0);
        fill(Bvec[1], -0.2, 0.4, 2.0);
        fill(Bvec[2], 0.5, -0.6, 3.0);
    }

    void applyAndCheckSpaced(BoundaryLocation loc, Box<std::uint32_t, 2> const& ghostBox)
    {
        auto B = Bvec.super();
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D> bc;
        bc.apply(B, loc, ghostBox, layout, 0.0);

        auto& Bx        = Bvec[0];
        auto& By        = Bvec[1];
        double const dx = layout.meshSize()[0];
        double const dy = layout.meshSize()[1];
        for (auto const& cell : ghostBox)
        {
            double const div
                = (Bx(cell.shift(0, 1)) - Bx(cell)) / dx + (By(cell.shift(1, 1)) - By(cell)) / dy;
            EXPECT_NEAR(div, 0.0, 1e-12)
                << "cell=(" << cell[0] << "," << cell[1] << ") at " << static_cast<int>(loc);
        }
    }
};

TEST_F(DivFreeTransverseNeumannBC2DAnisotropic, XBoundariesKillDivergenceOnAnisotropicMesh)
{
    applyAndCheckSpaced(BoundaryLocation::XLower, mhd2DXLowerGhostBox());
    applyAndCheckSpaced(BoundaryLocation::XUpper, mhd2DXUpperGhostBox());
}

TEST_F(DivFreeTransverseNeumannBC2DAnisotropic, YBoundariesKillDivergenceOnAnisotropicMesh)
{
    applyAndCheckSpaced(BoundaryLocation::YLower, mhd2DYLowerGhostBox());
    applyAndCheckSpaced(BoundaryLocation::YUpper, mhd2DYUpperGhostBox());
}


struct DivFreeTransverseNeumannBC3D : testing::Test
{
    GridLayoutMHD3D layout{
        {0.1, 0.1, 0.1}, {nCellsMHDX3D, nCellsMHDY3D, nCellsMHDZ3D}, {0.0, 0.0, 0.0}};


    UsableVecFieldMHD<3> Bvec{"bc_test_B", layout, MHDQuantity::Vector::B};


    DivFreeTransverseNeumannBC3D()
    {
        auto fill = [&](auto& f, double a, double b, double c, double d) {
            auto const shape = f.shape();
            for (std::uint32_t i = 0; i < shape[0]; ++i)
                for (std::uint32_t j = 0; j < shape[1]; ++j)
                    for (std::uint32_t k = 0; k < shape[2]; ++k)
                        f(i, j, k) = a * static_cast<double>(i) + b * static_cast<double>(j)
                                     + c * static_cast<double>(k) + d;
        };
        fill(Bvec[0], 0.3, 0.7, -0.2, 1.0);
        fill(Bvec[1], -0.2, 0.4, 0.5, 2.0);
        fill(Bvec[2], 0.5, -0.6, 0.3, 3.0);
    }

    void applyAndCheck(BoundaryLocation loc, Box<std::uint32_t, 3> const& ghostBox)
    {
        Direction const direction = getDirection(loc);
        Side const side           = getSide(loc);
        std::size_t const iNormal = static_cast<std::size_t>(direction);

        auto B = Bvec.super();
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D> bc;
        bc.apply(B, loc, ghostBox, layout, 0.0);

        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            if (comp == iNormal)
                continue;
            auto& bField         = Bvec[comp];
            auto const qty       = MHDQuantity::componentsQuantities(MHDQuantity::Vector::B)[comp];
            auto const centering = layout.centering(qty)[iNormal];
            for (auto const& index : layout.toFieldBox(ghostBox, qty))
            {
                auto const mirror = layout.boundaryMirrored(direction, side, centering, index);
                EXPECT_DOUBLE_EQ(bField(index), bField(mirror)) << "comp=" << comp;
            }
        }

        auto& Bx        = Bvec[0];
        auto& By        = Bvec[1];
        auto& Bz        = Bvec[2];
        double const dx = layout.meshSize()[0];
        double const dy = layout.meshSize()[1];
        double const dz = layout.meshSize()[2];
        for (auto const& cell : ghostBox)
        {
            double const div = (Bx(cell.shift(0, 1)) - Bx(cell)) / dx
                               + (By(cell.shift(1, 1)) - By(cell)) / dy
                               + (Bz(cell.shift(2, 1)) - Bz(cell)) / dz;
            EXPECT_NEAR(div, 0.0, 1e-12) << "at " << static_cast<int>(loc);
        }
    }
};

TEST_F(DivFreeTransverseNeumannBC3D, XBoundariesKillDivergence)
{
    applyAndCheck(BoundaryLocation::XLower, mhd3DXLowerGhostBox());
    applyAndCheck(BoundaryLocation::XUpper, mhd3DXUpperGhostBox());
}

TEST_F(DivFreeTransverseNeumannBC3D, YBoundariesKillDivergence)
{
    applyAndCheck(BoundaryLocation::YLower, mhd3DYLowerGhostBox());
    applyAndCheck(BoundaryLocation::YUpper, mhd3DYUpperGhostBox());
}

TEST_F(DivFreeTransverseNeumannBC3D, ZBoundariesKillDivergence)
{
    applyAndCheck(BoundaryLocation::ZLower, mhd3DZLowerGhostBox());
    applyAndCheck(BoundaryLocation::ZUpper, mhd3DZUpperGhostBox());
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
