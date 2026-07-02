#include "gtest/gtest.h"

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/field_adaptive_outflow_pressure_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_total_energy_from_pressure_boundary_condition.hpp"
#include "core/numerics/thermo/ideal_gas_thermo.hpp"
#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

using namespace PHARE::core;


// ════════════════════════════════════════════════════════════════════════════
// Adaptive outflow = FieldTotalEnergyFromPressureBoundaryCondition with a
// FieldAdaptiveOutflowPressureBoundaryCondition as its P_bc.
//
// Per tangential column, the pressure ghost is chosen from the first-interior
// regime (u_n = sign_n · V_n vs the boundary-normal fast magnetosonic speed c_f):
//   * super-magnetofast (u_n ≥ c_f): pressure is zero-gradient → reconstructed
//     energy equals the interior energy (etot_val).
//   * sub-fast (u_n < c_f): pressure is Dirichlet(P_fixed) → energy reconstructed
//     from the fixed ghost pressure (etot_ghost), same as a fixed-pressure outlet.
//
// A single large +x velocity makes XUpper super-fast (u_n = +vx) while XLower and
// the transverse sides (u_n = ±vy/±vz with vy=vz=0) stay sub-fast — exercising
// both branches in one fixture.
//
// Non-corner ghost-box helpers exclude cells overwritten by adjacent-direction
// BCs (corner values are order-dependent and not read by physical stencils).
// ════════════════════════════════════════════════════════════════════════════

static Box<std::uint32_t, 2> ao2DXLowerNonCornerGhostBox()
{
    return {{0u, mhdGhostWidth}, {mhdGhostWidth - 1u, mhdGhostWidth + nCellsMHDY2D - 1u}};
}
static Box<std::uint32_t, 2> ao2DXUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth + nCellsMHDX2D, mhdGhostWidth},
            {2u * mhdGhostWidth + nCellsMHDX2D - 1u, mhdGhostWidth + nCellsMHDY2D - 1u}};
}
static Box<std::uint32_t, 2> ao2DYLowerNonCornerGhostBox()
{
    return {{mhdGhostWidth, 0u}, {mhdGhostWidth + nCellsMHDX2D - 1u, mhdGhostWidth - 1u}};
}
static Box<std::uint32_t, 2> ao2DYUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth, mhdGhostWidth + nCellsMHDY2D},
            {mhdGhostWidth + nCellsMHDX2D - 1u, 2u * mhdGhostWidth + nCellsMHDY2D - 1u}};
}

static Box<std::uint32_t, 3> ao3DXLowerNonCornerGhostBox()
{
    return {{0u, mhdGhostWidth, mhdGhostWidth},
            {mhdGhostWidth - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> ao3DXUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth + nCellsMHDX3D, mhdGhostWidth, mhdGhostWidth},
            {2u * mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> ao3DYLowerNonCornerGhostBox()
{
    return {{mhdGhostWidth, 0u, mhdGhostWidth},
            {mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> ao3DZLowerNonCornerGhostBox()
{
    return {{mhdGhostWidth, mhdGhostWidth, 0u},
            {mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth - 1u}};
}


// Shared physical state: large +x velocity → XUpper super-fast, all other sides sub-fast.
struct AdaptiveOutflowState
{
    static constexpr double gamma   = 5.0 / 3.0;
    static constexpr double rho_val = 2.0;
    static constexpr double vx_val  = 10.0; // ≫ c_f → super-fast in +x
    static constexpr double vy_val  = 0.0;
    static constexpr double vz_val  = 0.0;
    static constexpr double Bx_val  = 0.5;
    static constexpr double By_val  = 0.5;
    static constexpr double Bz_val  = 0.5;
    static constexpr double P_val   = 1.0;
    static constexpr double P_fixed = 0.5;

    static constexpr double u_specific_interior = P_val / (rho_val * (gamma - 1.0));
    static constexpr double etot_val
        = rho_val * u_specific_interior
          + 0.5 * rho_val * (vx_val * vx_val + vy_val * vy_val + vz_val * vz_val)
          + 0.5 * (Bx_val * Bx_val + By_val * By_val + Bz_val * Bz_val);

    static constexpr double P_ghost          = 2.0 * P_fixed - P_val;
    static constexpr double u_specific_ghost = P_ghost / (rho_val * (gamma - 1.0));
    static constexpr double etot_ghost
        = rho_val * u_specific_ghost
          + 0.5 * rho_val * (vx_val * vx_val + vy_val * vy_val + vz_val * vz_val)
          + 0.5 * (Bx_val * Bx_val + By_val * By_val + Bz_val * Bz_val);
};


// ─── 1D ──────────────────────────────────────────────────────────────────────

struct AdaptiveOutflowBC1D : testing::Test, AdaptiveOutflowState
{
    GridLayoutMHD1D layout{{0.1}, {nCellsMHD}, {0.0}};

    GridMHD1D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD1D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD1D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<1> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<1> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<1> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<1>& rhoField{*(&rhoGrid)};
    FieldMHD<1>& EtotField{*(&EtotGrid)};

    AdaptiveOutflowBC1D()
    {
        auto fill_all = [&](auto& f, double val) {
            for (std::uint32_t i = 0; i < f.shape()[0]; ++i)
                f(i) = val;
        };
        fill_all(rhoField, rho_val);
        fill_all(PGrid, P_val);
        fill_all(EtotField, etot_val);
        fill_all(rhoV[0], rho_val * vx_val);
        fill_all(rhoV[1], rho_val * vy_val);
        fill_all(rhoV[2], rho_val * vz_val);
        fill_all(Bvec[0], Bx_val);
        fill_all(Bvec[1], By_val);
        fill_all(Bvec[2], Bz_val);
    }

    auto makeBC()
    {
        auto rho_bc
            = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
        auto rhoV_bc
            = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
        auto B_bc = std::make_shared<
            FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
        auto thermo = std::make_shared<IdealGasThermo>(gamma);
        auto P_bc   = std::make_shared<
            FieldAdaptiveOutflowPressureBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>(P_fixed,
                                                                                        thermo);
        return std::make_pair(
            FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>{
                rho_bc, rhoV_bc, B_bc, P_bc, thermo},
            thermo);
    }
};

TEST_F(AdaptiveOutflowBC1D, SuperFastUpperIsZeroGradient_SubFastLowerIsFixedPressure)
{
    auto [bc, thermo] = makeBC();

    bc.apply(EtotField, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhdUpperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto etotQty     = MHDQuantity::Scalar::Etot;
    std::uint32_t ps = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe = layout.physicalEndIndex(etotQty, Direction::X);

    for (std::uint32_t g = 0; g < mhdGhostWidth; ++g)
    {
        EXPECT_NEAR(EtotField(ps - 1u - g), etot_ghost, 1e-12) << "XLower (sub-fast) g=" << g;
        EXPECT_NEAR(EtotField(pe + 1u + g), etot_val, 1e-12) << "XUpper (super-fast) g=" << g;
    }
}


// ─── 2D ──────────────────────────────────────────────────────────────────────

struct AdaptiveOutflowBC2D : testing::Test, AdaptiveOutflowState
{
    GridLayoutMHD2D layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};

    GridMHD2D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD2D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD2D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<2> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<2> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<2> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<2>& EtotField{*(&EtotGrid)};

    AdaptiveOutflowBC2D()
    {
        auto fill_all = [&](auto& f, double val) {
            auto [nx, ny] = f.shape();
            for (std::uint32_t ix = 0; ix < nx; ++ix)
                for (std::uint32_t iy = 0; iy < ny; ++iy)
                    f(ix, iy) = val;
        };
        fill_all(rhoGrid, rho_val);
        fill_all(PGrid, P_val);
        fill_all(EtotField, etot_val);
        fill_all(rhoV[0], rho_val * vx_val);
        fill_all(rhoV[1], rho_val * vy_val);
        fill_all(rhoV[2], rho_val * vz_val);
        fill_all(Bvec[0], Bx_val);
        fill_all(Bvec[1], By_val);
        fill_all(Bvec[2], Bz_val);
    }

    auto makeBC()
    {
        auto rho_bc
            = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
        auto rhoV_bc
            = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
        auto B_bc = std::make_shared<
            FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
        auto thermo = std::make_shared<IdealGasThermo>(gamma);
        auto P_bc   = std::make_shared<
            FieldAdaptiveOutflowPressureBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>(P_fixed,
                                                                                        thermo);
        return std::make_pair(
            FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>{
                rho_bc, rhoV_bc, B_bc, P_bc, thermo},
            thermo);
    }
};

TEST_F(AdaptiveOutflowBC2D, SuperFastUpperIsZeroGradient_OtherSidesFixedPressure)
{
    auto [bc, thermo] = makeBC();

    bc.apply(EtotField, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd2DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd2DYUpperGhostBox(), layout, makeCtx(acc, 0.0));

    for (auto const& idx : ao2DXUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "XUpper super-fast (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : ao2DXLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "XLower sub-fast (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : ao2DYLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YLower sub-fast (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : ao2DYUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YUpper sub-fast (" << idx[0] << "," << idx[1] << ")";
}


// ─── 3D ──────────────────────────────────────────────────────────────────────

struct AdaptiveOutflowBC3D : testing::Test, AdaptiveOutflowState
{
    GridLayoutMHD3D layout{
        {0.1, 0.1, 0.1}, {nCellsMHDX3D, nCellsMHDY3D, nCellsMHDZ3D}, {0.0, 0.0, 0.0}};

    GridMHD3D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD3D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD3D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<3> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<3> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<3> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<3>& EtotField{*(&EtotGrid)};

    AdaptiveOutflowBC3D()
    {
        auto fill_all = [&](auto& f, double val) {
            auto [nx, ny, nz] = f.shape();
            for (std::uint32_t ix = 0; ix < nx; ++ix)
                for (std::uint32_t iy = 0; iy < ny; ++iy)
                    for (std::uint32_t iz = 0; iz < nz; ++iz)
                        f(ix, iy, iz) = val;
        };
        fill_all(rhoGrid, rho_val);
        fill_all(PGrid, P_val);
        fill_all(EtotField, etot_val);
        fill_all(rhoV[0], rho_val * vx_val);
        fill_all(rhoV[1], rho_val * vy_val);
        fill_all(rhoV[2], rho_val * vz_val);
        fill_all(Bvec[0], Bx_val);
        fill_all(Bvec[1], By_val);
        fill_all(Bvec[2], Bz_val);
    }

    auto makeBC()
    {
        auto rho_bc
            = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
        auto rhoV_bc
            = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
        auto B_bc = std::make_shared<
            FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
        auto thermo = std::make_shared<IdealGasThermo>(gamma);
        auto P_bc   = std::make_shared<
            FieldAdaptiveOutflowPressureBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>(P_fixed,
                                                                                        thermo);
        return std::make_pair(
            FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>{
                rho_bc, rhoV_bc, B_bc, P_bc, thermo},
            thermo);
    }
};

TEST_F(AdaptiveOutflowBC3D, SuperFastUpperIsZeroGradient_OtherSidesFixedPressure)
{
    auto [bc, thermo] = makeBC();

    bc.apply(EtotField, BoundaryLocation::XLower, mhd3DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd3DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd3DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd3DYUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZLower, mhd3DZLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZUpper, mhd3DZUpperGhostBox(), layout, makeCtx(acc, 0.0));

    for (auto const& idx : ao3DXUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "XUpper super-fast (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : ao3DXLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "XLower sub-fast (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : ao3DYLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YLower sub-fast (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : ao3DZLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "ZLower sub-fast (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
