#include "gtest/gtest.h"

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_neumann_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

using namespace PHARE::core;


// ════════════════════════════════════════════════════════════════════════════
// FieldDivergenceFreeTransverseNeumannBoundaryCondition
//
// Transverse components mirror the first interior value (zero normal gradient),
// the normal component is set so that the discrete divergence of B vanishes in
// every ghost cell.
// ════════════════════════════════════════════════════════════════════════════

struct DivFreeTransverseNeumannBC2D : testing::Test
{
    GridLayoutMHD2D layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};

    GridMHD2D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD2D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD2D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<2> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<2> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<2> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    DivFreeTransverseNeumannBC2D()
    {
        // analytic node values everywhere (interior *and* ghosts), deliberately not
        // divergence-free so the BC has real work to do on the normal component
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
                = (Bx(cell.neighbor(0, 1)) - Bx(cell)) + (By(cell.neighbor(1, 1)) - By(cell));
            EXPECT_NEAR(div, 0.0, 1e-12)
                << "cell=(" << cell[0] << "," << cell[1] << ") at " << static_cast<int>(loc);
        }
    }

    void applyAndCheck(BoundaryLocation loc, Box<std::uint32_t, 2> const& ghostBox)
    {
        auto B = Bvec.super();
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D> bc;
        bc.apply(B, loc, ghostBox, layout, makeCtx(acc));

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


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
