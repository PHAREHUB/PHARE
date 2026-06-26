#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "core/numerics/MHD_equations/MHD_equations.hpp"
#include "core/numerics/riemann_solvers/rusanov.hpp"

#include "gtest/gtest.h"

using namespace PHARE::core;

// Well-balancing of the MHD momentum flux for a static background field B0.
//
// The flux is built from the total field B = B0 + B1. With the static-stress
// subtraction in MHDEquations::compute, a rest state (B1 = 0, V = 0, uniform P)
// must produce a momentum flux that depends ONLY on the thermal pressure, i.e.
// (P0, 0, 0) on every face regardless of the reconstructed B0 left/right values.
// The discrete divergence across faces is then exactly zero -> machine-precision
// steady state. Without the subtraction the magnetic stress of B0 leaks into the
// flux and these expectations fail.
namespace
{
constexpr double gamma = 5.0 / 3.0;
constexpr double rho0  = 1.5;
constexpr double P0    = 2.0;

// A rest PerIndex: perturbation B1 = 0 (PerIndex stores B1), background B0, V = 0.
auto restState(PerIndexVector<double> const& B0)
{
    return PerIndex<double>{rho0, PerIndexVector<double>{0.0, 0.0, 0.0},
                            PerIndexVector<double>{0.0, 0.0, 0.0}, P0, B0};
}

template<auto direction>
auto solveFace(PerIndexVector<double> const& B0L, PerIndexVector<double> const& B0R)
{
    MHDEquations<false, false, false> equations{gamma, 0.0, 0.0};
    Rusanov<false> riemann{gamma};

    auto uL = restState(B0L);
    auto uR = restState(B0R);

    auto const fL = equations.template compute<direction>(uL);
    auto const fR = equations.template compute<direction>(uR);

    // solve mutates uL/uR to conservative form internally; fluxes computed above.
    return riemann.template solve<direction>(uL, uR, fL, fR);
}
} // namespace

TEST(MHDWellBalanced, momentumFluxIsBackgroundIndependentX)
{
    auto const flux = solveFace<Direction::X>({0.3, 0.7, -0.2}, {0.5, -0.4, 0.6});

    EXPECT_DOUBLE_EQ(flux.rhoV().x, P0);
    EXPECT_DOUBLE_EQ(flux.rhoV().y, 0.0);
    EXPECT_DOUBLE_EQ(flux.rhoV().z, 0.0);
}

TEST(MHDWellBalanced, momentumFluxIsBackgroundIndependentY)
{
    auto const flux = solveFace<Direction::Y>({0.3, 0.7, -0.2}, {0.5, -0.4, 0.6});

    EXPECT_DOUBLE_EQ(flux.rhoV().x, 0.0);
    EXPECT_DOUBLE_EQ(flux.rhoV().y, P0);
    EXPECT_DOUBLE_EQ(flux.rhoV().z, 0.0);
}

// Two adjacent faces carrying different reconstructed B0 must yield identical
// momentum flux -> the finite-volume divergence (and thus the rhoV update) is
// exactly zero at rest.
TEST(MHDWellBalanced, restStateHasZeroMomentumDivergenceX)
{
    auto const faceLeft  = solveFace<Direction::X>({0.3, 0.7, -0.2}, {0.5, -0.4, 0.6});
    auto const faceRight = solveFace<Direction::X>({0.5, -0.4, 0.6}, {1.1, 0.2, 0.9});

    EXPECT_DOUBLE_EQ(faceRight.rhoV().x - faceLeft.rhoV().x, 0.0);
    EXPECT_DOUBLE_EQ(faceRight.rhoV().y - faceLeft.rhoV().y, 0.0);
    EXPECT_DOUBLE_EQ(faceRight.rhoV().z - faceLeft.rhoV().z, 0.0);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
