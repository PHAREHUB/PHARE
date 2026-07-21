#include "gtest/gtest.h"

#include "core/numerics/boundary_condition/field_divergence_free_transverse_dirichlet_boundary_condition.hpp"
#include "tests/core/numerics/boundary_condition/hybrid_bc_test_fixtures.hpp"

using namespace PHARE::core;


// ─── 1D VecField ─────────────────────────────────────────────────────────────

TEST_F(VecFieldBC1D, DivergenceFreeTransverseDirichletAtXBoundaries)
{
    std::array values{123.0, 7.0, 11.0};
    FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecField1D, GridLayout1D> bc{values};
    bc.apply(B, BoundaryLocation::XLower, lowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::XUpper, upperGhostCellBox(), layout, makeCtx(acc, 0.0));

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


// ─── 2D VecField ─────────────────────────────────────────────────────────────

TEST_F(VecFieldBC2D, DivergenceFreeTransverseDirichletAtXBoundaries)
{
    std::array values{123.0, 7.0, 11.0};
    FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecField2D, GridLayout2D> bc{values};
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, makeCtx(acc, 0.0));

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
    FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecField2D, GridLayout2D> bc{values};
    bc.apply(B, BoundaryLocation::YLower, yLowerGhostCellBox2D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::YUpper, yUpperGhostCellBox2D(), layout, makeCtx(acc, 0.0));

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
    FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecField2D, GridLayout2D> bc{
        std::array{123.0, 0.0, 11.0}};
    bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, makeCtx(acc, 0.0));

    auto& Bx = B[0];
    auto& By = B[1];
    for (auto const& index : xLowerGhostCellBox2D())
    {
        EXPECT_DOUBLE_EQ(Bx(index.template neighbor<0, 1>()) - Bx(index)
                             + By(index.template neighbor<1, 1>()) - By(index),
                         0.0)
            << "lower divergence at (" << index[0] << ", " << index[1] << ")";
    }
    for (auto const& index : xUpperGhostCellBox2D())
    {
        EXPECT_DOUBLE_EQ(Bx(index.template neighbor<0, 1>()) - Bx(index)
                             + By(index.template neighbor<1, 1>()) - By(index),
                         0.0)
            << "upper divergence at (" << index[0] << ", " << index[1] << ")";
    }
}


// ─── 3D VecField ─────────────────────────────────────────────────────────────

TEST_F(VecFieldBC3D, DivergenceFreeTransverseDirichletAtZBoundaries)
{
    std::array values{3.0, 7.0, 123.0};
    FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecField3D, GridLayout3D> bc{values};
    bc.apply(B, BoundaryLocation::ZLower, zLowerGhostCellBox3D(), layout, makeCtx(acc, 0.0));
    bc.apply(B, BoundaryLocation::ZUpper, zUpperGhostCellBox3D(), layout, makeCtx(acc, 0.0));

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


// Time-varying prescribed tangential values: each tangential component is driven by a
// space-time function evaluated at ctx.time. With spatially-uniform functions the tangential
// ghosts extrapolate to 2*f(t) - interior and the normal component stays divergence-free.
// Running at two times proves ctx.time threads through the per-component scalar Dirichlets.
TEST_F(VecFieldBC2D, DivergenceFreeTransverseDirichletTimeVaryingAtXBoundaries)
{
    auto uniformFn = [](double base) {
        return PHARE::initializer::SpaceTimeFunction<2>{
            [base](std::vector<double> const& x, std::vector<double> const& /*y*/, double t) {
                std::vector<double> out(x.size(), base + t);
                return std::shared_ptr<Span<double>>{
                    std::make_shared<VectorSpan<double>>(std::move(out))};
            }};
    };

    // component 0 (Bx) is the X-normal component: its function is never applied (normal comes
    // from divB=0), only the tangential components 1 (By) and 2 (Bz) are Dirichlet-driven.
    std::array<double, 3> const base{0.0, 7.0, 11.0};
    std::array<PHARE::initializer::SpaceTimeFunction<2>, 3> fns{uniformFn(base[0]),
                                                               uniformFn(base[1]),
                                                               uniformFn(base[2])};
    FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecField2D, GridLayout2D> bc{fns};

    auto checkAtTime = [&](double const t) {
        // reset interior to the constant reference value on every pass
        for (std::size_t comp = 0; comp < 3; ++comp)
        {
            auto& f          = B[comp];
            auto qty         = HybridQuantity::componentsQuantities(vecQty)[comp];
            std::uint32_t sx = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t ex = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t sy = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t ey = layout.physicalEndIndex(qty, Direction::Y);
            for (std::uint32_t ix = sx; ix <= ex; ++ix)
                for (std::uint32_t iy = sy; iy <= ey; ++iy)
                    f(ix, iy) = interiorValue;
        }

        bc.apply(B, BoundaryLocation::XLower, xLowerGhostCellBox2D(), layout, makeCtx(acc, t));
        bc.apply(B, BoundaryLocation::XUpper, xUpperGhostCellBox2D(), layout, makeCtx(acc, t));

        for (std::size_t comp : {1u, 2u})
        {
            auto& f              = B[comp];
            auto qty             = HybridQuantity::componentsQuantities(vecQty)[comp];
            std::uint32_t psx    = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t pex    = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t sy     = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t ey     = layout.physicalEndIndex(qty, Direction::Y);
            double expectedGhost = 2.0 * (base[comp] + t) - interiorValue;
            for (std::uint32_t iy = sy; iy <= ey; ++iy)
            {
                EXPECT_DOUBLE_EQ(f(psx - 1, iy), expectedGhost)
                    << "comp=" << comp << " lower near ghost iy=" << iy << " t=" << t;
                EXPECT_DOUBLE_EQ(f(pex + 1, iy), expectedGhost)
                    << "comp=" << comp << " upper near ghost iy=" << iy << " t=" << t;
            }
        }
    };

    checkAtTime(0.0);
    checkAtTime(2.5);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
