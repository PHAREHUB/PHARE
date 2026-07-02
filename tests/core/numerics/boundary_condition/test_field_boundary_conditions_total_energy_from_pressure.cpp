#include "gtest/gtest.h"

#include "core/boundary/boundary_defs.hpp"
#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_total_energy_from_pressure_boundary_condition.hpp"
#include "core/numerics/thermo/ideal_gas_thermo.hpp"
#include "tests/core/numerics/boundary_condition/mhd_bc_test_fixtures.hpp"

using namespace PHARE::core;


// ════════════════════════════════════════════════════════════════════════════
// FieldTotalEnergyFromPressureBoundaryCondition — all-Neumann sub-BCs
//
// When ρ, ρv, B, P all use Neumann, the ghost thermodynamic state mirrors the
// interior → ghost Etot equals interior Etot.
// ════════════════════════════════════════════════════════════════════════════

// ─── 1D ─────────────────────────────────────────────────────────────────────

struct EtotFromPressureBC1D : testing::Test
{
    static constexpr double gamma    = 5.0 / 3.0;
    static constexpr double rho_val  = 2.0;
    static constexpr double vx_val   = 1.0;
    static constexpr double vy_val   = 1.0;
    static constexpr double vz_val   = 1.0;
    static constexpr double Bx_val   = 0.5;
    static constexpr double By_val   = 0.5;
    static constexpr double Bz_val   = 0.5;
    static constexpr double P_val    = 1.0;
    static constexpr double sentinel = -999.0;

    static constexpr double u_specific = P_val / (rho_val * (gamma - 1.0));
    static constexpr double etot_val
        = rho_val * u_specific
          + 0.5 * rho_val * (vx_val * vx_val + vy_val * vy_val + vz_val * vz_val)
          + 0.5 * (Bx_val * Bx_val + By_val * By_val + Bz_val * Bz_val);

    GridLayoutMHD1D layout{{0.1}, {nCellsMHD}, {0.0}};

    GridMHD1D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD1D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD1D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<1> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<1> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<1> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<1>& rhoField{*(&rhoGrid)};
    FieldMHD<1>& PField{*(&PGrid)};
    FieldMHD<1>& EtotField{*(&EtotGrid)};

    EtotFromPressureBC1D()
    {
        auto fill_scalar = [&](auto& f, MHDQuantity::Scalar qty, double interior_val) {
            for (std::uint32_t i = 0; i < f.shape()[0]; ++i)
                f(i) = sentinel;
            std::uint32_t ps = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t pe = layout.physicalEndIndex(qty, Direction::X);
            for (std::uint32_t i = ps; i <= pe; ++i)
                f(i) = interior_val;
        };

        fill_scalar(rhoField, MHDQuantity::Scalar::rho, rho_val);
        fill_scalar(PField, MHDQuantity::Scalar::P, P_val);
        fill_scalar(EtotField, MHDQuantity::Scalar::Etot, etot_val);

        fill_scalar(rhoV[0], MHDQuantity::Scalar::rhoVx, rho_val * vx_val);
        fill_scalar(rhoV[1], MHDQuantity::Scalar::rhoVy, rho_val * vy_val);
        fill_scalar(rhoV[2], MHDQuantity::Scalar::rhoVz, rho_val * vz_val);

        fill_scalar(Bvec[0], MHDQuantity::Scalar::Bx, Bx_val);
        fill_scalar(Bvec[1], MHDQuantity::Scalar::By, By_val);
        fill_scalar(Bvec[2], MHDQuantity::Scalar::Bz, Bz_val);
    }
};

TEST_F(EtotFromPressureBC1D, NeumannSubBCsGhostEtotEqualsInteriorEtot)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
    auto P_bc   = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto B_bc = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhdUpperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto etotQty     = MHDQuantity::Scalar::Etot;
    std::uint32_t ps = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe = layout.physicalEndIndex(etotQty, Direction::X);

    for (std::uint32_t g = 0; g < mhdGhostWidth; ++g)
    {
        EXPECT_NEAR(EtotField(ps - 1u - g), etot_val, 1e-12) << "lower ghost g=" << g;
        EXPECT_NEAR(EtotField(pe + 1u + g), etot_val, 1e-12) << "upper ghost g=" << g;
    }
}

TEST_F(EtotFromPressureBC1D, InteriorEtotUnchangedAfterBC)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
    auto P_bc   = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto B_bc = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhdUpperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto etotQty     = MHDQuantity::Scalar::Etot;
    std::uint32_t ps = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe = layout.physicalEndIndex(etotQty, Direction::X);

    for (std::uint32_t i = ps; i <= pe; ++i)
        EXPECT_DOUBLE_EQ(EtotField(i), etot_val) << "interior index i=" << i;
}


// ─── 2D ─────────────────────────────────────────────────────────────────────

struct EtotFromPressureBC2D : testing::Test
{
    static constexpr double gamma    = 5.0 / 3.0;
    static constexpr double rho_val  = 2.0;
    static constexpr double vx_val   = 1.0;
    static constexpr double vy_val   = 1.0;
    static constexpr double vz_val   = 1.0;
    static constexpr double Bx_val   = 0.5;
    static constexpr double By_val   = 0.5;
    static constexpr double Bz_val   = 0.5;
    static constexpr double P_val    = 1.0;
    static constexpr double sentinel = -999.0;

    static constexpr double u_specific = P_val / (rho_val * (gamma - 1.0));
    static constexpr double etot_val
        = rho_val * u_specific
          + 0.5 * rho_val * (vx_val * vx_val + vy_val * vy_val + vz_val * vz_val)
          + 0.5 * (Bx_val * Bx_val + By_val * By_val + Bz_val * Bz_val);

    GridLayoutMHD2D layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};

    GridMHD2D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD2D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD2D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<2> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<2> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<2> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<2>& rhoField{*(&rhoGrid)};
    FieldMHD<2>& PField{*(&PGrid)};
    FieldMHD<2>& EtotField{*(&EtotGrid)};

    EtotFromPressureBC2D()
    {
        auto fill_scalar = [&](auto& f, MHDQuantity::Scalar qty, double interior_val) {
            auto const sz = f.shape();
            for (std::uint32_t ix = 0; ix < sz[0]; ++ix)
                for (std::uint32_t iy = 0; iy < sz[1]; ++iy)
                    f(ix, iy) = sentinel;
            std::uint32_t ps_x = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t pe_x = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t ps_y = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t pe_y = layout.physicalEndIndex(qty, Direction::Y);
            for (std::uint32_t ix = ps_x; ix <= pe_x; ++ix)
                for (std::uint32_t iy = ps_y; iy <= pe_y; ++iy)
                    f(ix, iy) = interior_val;
        };

        fill_scalar(rhoField, MHDQuantity::Scalar::rho, rho_val);
        fill_scalar(PField, MHDQuantity::Scalar::P, P_val);
        fill_scalar(EtotField, MHDQuantity::Scalar::Etot, etot_val);

        fill_scalar(rhoV[0], MHDQuantity::Scalar::rhoVx, rho_val * vx_val);
        fill_scalar(rhoV[1], MHDQuantity::Scalar::rhoVy, rho_val * vy_val);
        fill_scalar(rhoV[2], MHDQuantity::Scalar::rhoVz, rho_val * vz_val);

        fill_scalar(Bvec[0], MHDQuantity::Scalar::Bx, Bx_val);
        fill_scalar(Bvec[1], MHDQuantity::Scalar::By, By_val);
        fill_scalar(Bvec[2], MHDQuantity::Scalar::Bz, Bz_val);
    }
};

TEST_F(EtotFromPressureBC2D, NeumannSubBCsGhostEtotEqualsInteriorEtot)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
    auto P_bc   = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto B_bc = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd2DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd2DYUpperGhostBox(), layout, makeCtx(acc, 0.0));

    for (auto const& idx : mhd2DXLowerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "XLower ghost (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : mhd2DXUpperGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "XUpper ghost (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : mhd2DYLowerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "YLower ghost (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : mhd2DYUpperGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "YUpper ghost (" << idx[0] << "," << idx[1] << ")";
}

TEST_F(EtotFromPressureBC2D, InteriorEtotUnchangedAfterBC)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
    auto P_bc   = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto B_bc = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd2DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd2DYUpperGhostBox(), layout, makeCtx(acc, 0.0));

    auto etotQty       = MHDQuantity::Scalar::Etot;
    std::uint32_t ps_x = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe_x = layout.physicalEndIndex(etotQty, Direction::X);
    std::uint32_t ps_y = layout.physicalStartIndex(etotQty, Direction::Y);
    std::uint32_t pe_y = layout.physicalEndIndex(etotQty, Direction::Y);

    for (std::uint32_t ix = ps_x; ix <= pe_x; ++ix)
        for (std::uint32_t iy = ps_y; iy <= pe_y; ++iy)
            EXPECT_DOUBLE_EQ(EtotField(ix, iy), etot_val)
                << "interior index (" << ix << "," << iy << ")";
}


// ─── 3D ─────────────────────────────────────────────────────────────────────

struct EtotFromPressureBC3D : testing::Test
{
    static constexpr double gamma    = 5.0 / 3.0;
    static constexpr double rho_val  = 2.0;
    static constexpr double vx_val   = 1.0;
    static constexpr double vy_val   = 1.0;
    static constexpr double vz_val   = 1.0;
    static constexpr double Bx_val   = 0.5;
    static constexpr double By_val   = 0.5;
    static constexpr double Bz_val   = 0.5;
    static constexpr double P_val    = 1.0;
    static constexpr double sentinel = -999.0;

    static constexpr double u_specific = P_val / (rho_val * (gamma - 1.0));
    static constexpr double etot_val
        = rho_val * u_specific
          + 0.5 * rho_val * (vx_val * vx_val + vy_val * vy_val + vz_val * vz_val)
          + 0.5 * (Bx_val * Bx_val + By_val * By_val + Bz_val * Bz_val);

    GridLayoutMHD3D layout{
        {0.1, 0.1, 0.1}, {nCellsMHDX3D, nCellsMHDY3D, nCellsMHDZ3D}, {0.0, 0.0, 0.0}};

    GridMHD3D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD3D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD3D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<3> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<3> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<3> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<3>& rhoField{*(&rhoGrid)};
    FieldMHD<3>& PField{*(&PGrid)};
    FieldMHD<3>& EtotField{*(&EtotGrid)};

    EtotFromPressureBC3D()
    {
        auto fill_scalar = [&](auto& f, MHDQuantity::Scalar qty, double interior_val) {
            auto const sz = f.shape();
            for (std::uint32_t ix = 0; ix < sz[0]; ++ix)
                for (std::uint32_t iy = 0; iy < sz[1]; ++iy)
                    for (std::uint32_t iz = 0; iz < sz[2]; ++iz)
                        f(ix, iy, iz) = sentinel;
            std::uint32_t ps_x = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t pe_x = layout.physicalEndIndex(qty, Direction::X);
            std::uint32_t ps_y = layout.physicalStartIndex(qty, Direction::Y);
            std::uint32_t pe_y = layout.physicalEndIndex(qty, Direction::Y);
            std::uint32_t ps_z = layout.physicalStartIndex(qty, Direction::Z);
            std::uint32_t pe_z = layout.physicalEndIndex(qty, Direction::Z);
            for (std::uint32_t ix = ps_x; ix <= pe_x; ++ix)
                for (std::uint32_t iy = ps_y; iy <= pe_y; ++iy)
                    for (std::uint32_t iz = ps_z; iz <= pe_z; ++iz)
                        f(ix, iy, iz) = interior_val;
        };

        fill_scalar(rhoField, MHDQuantity::Scalar::rho, rho_val);
        fill_scalar(PField, MHDQuantity::Scalar::P, P_val);
        fill_scalar(EtotField, MHDQuantity::Scalar::Etot, etot_val);

        fill_scalar(rhoV[0], MHDQuantity::Scalar::rhoVx, rho_val * vx_val);
        fill_scalar(rhoV[1], MHDQuantity::Scalar::rhoVy, rho_val * vy_val);
        fill_scalar(rhoV[2], MHDQuantity::Scalar::rhoVz, rho_val * vz_val);

        fill_scalar(Bvec[0], MHDQuantity::Scalar::Bx, Bx_val);
        fill_scalar(Bvec[1], MHDQuantity::Scalar::By, By_val);
        fill_scalar(Bvec[2], MHDQuantity::Scalar::Bz, Bz_val);
    }
};

TEST_F(EtotFromPressureBC3D, NeumannSubBCsGhostEtotEqualsInteriorEtot)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
    auto P_bc   = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto B_bc = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd3DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd3DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd3DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd3DYUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZLower, mhd3DZLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZUpper, mhd3DZUpperGhostBox(), layout, makeCtx(acc, 0.0));

    for (auto const& idx : mhd3DXLowerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "XLower ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : mhd3DXUpperGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "XUpper ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : mhd3DYLowerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "YLower ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : mhd3DYUpperGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "YUpper ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : mhd3DZLowerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "ZLower ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : mhd3DZUpperGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_val, 1e-12)
            << "ZUpper ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
}

TEST_F(EtotFromPressureBC3D, InteriorEtotUnchangedAfterBC)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
    auto P_bc   = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto B_bc = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd3DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd3DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd3DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd3DYUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZLower, mhd3DZLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZUpper, mhd3DZUpperGhostBox(), layout, makeCtx(acc, 0.0));

    auto etotQty       = MHDQuantity::Scalar::Etot;
    std::uint32_t ps_x = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe_x = layout.physicalEndIndex(etotQty, Direction::X);
    std::uint32_t ps_y = layout.physicalStartIndex(etotQty, Direction::Y);
    std::uint32_t pe_y = layout.physicalEndIndex(etotQty, Direction::Y);
    std::uint32_t ps_z = layout.physicalStartIndex(etotQty, Direction::Z);
    std::uint32_t pe_z = layout.physicalEndIndex(etotQty, Direction::Z);

    for (std::uint32_t ix = ps_x; ix <= pe_x; ++ix)
        for (std::uint32_t iy = ps_y; iy <= pe_y; ++iy)
            for (std::uint32_t iz = ps_z; iz <= pe_z; ++iz)
                EXPECT_DOUBLE_EQ(EtotField(ix, iy, iz), etot_val)
                    << "interior index (" << ix << "," << iy << "," << iz << ")";
}


// ════════════════════════════════════════════════════════════════════════════
// FieldTotalEnergyFromPressureBoundaryCondition — FixedPressureOutflow variant
//
// ρ, ρv use Neumann; B uses DivergenceFreeTransverseNeumann; P uses Dirichlet(P_fixed).
// Ghost Etot is reconstructed from the fixed ghost pressure + Neumann-extrapolated
// density, momentum, and magnetic field.
//
// Non-corner ghost-box helpers exclude cells overwritten by adjacent-direction
// BCs (corner values are order-dependent and not read by physical stencils).
// ════════════════════════════════════════════════════════════════════════════

static Box<std::uint32_t, 2> fp2DXLowerNonCornerGhostBox()
{
    return {{0u, mhdGhostWidth},
            {mhdGhostWidth - 1u, mhdGhostWidth + nCellsMHDY2D - 1u}};
}
static Box<std::uint32_t, 2> fp2DXUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth + nCellsMHDX2D, mhdGhostWidth},
            {2u * mhdGhostWidth + nCellsMHDX2D - 1u, mhdGhostWidth + nCellsMHDY2D - 1u}};
}
static Box<std::uint32_t, 2> fp2DYLowerNonCornerGhostBox()
{
    return {{mhdGhostWidth, 0u},
            {mhdGhostWidth + nCellsMHDX2D - 1u, mhdGhostWidth - 1u}};
}
static Box<std::uint32_t, 2> fp2DYUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth, mhdGhostWidth + nCellsMHDY2D},
            {mhdGhostWidth + nCellsMHDX2D - 1u, 2u * mhdGhostWidth + nCellsMHDY2D - 1u}};
}

static Box<std::uint32_t, 3> fp3DXLowerNonCornerGhostBox()
{
    return {{0u, mhdGhostWidth, mhdGhostWidth},
            {mhdGhostWidth - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> fp3DXUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth + nCellsMHDX3D, mhdGhostWidth, mhdGhostWidth},
            {2u * mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> fp3DYLowerNonCornerGhostBox()
{
    return {{mhdGhostWidth, 0u, mhdGhostWidth},
            {mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> fp3DYUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth, mhdGhostWidth + nCellsMHDY3D, mhdGhostWidth},
            {mhdGhostWidth + nCellsMHDX3D - 1u, 2u * mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth + nCellsMHDZ3D - 1u}};
}
static Box<std::uint32_t, 3> fp3DZLowerNonCornerGhostBox()
{
    return {{mhdGhostWidth, mhdGhostWidth, 0u},
            {mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             mhdGhostWidth - 1u}};
}
static Box<std::uint32_t, 3> fp3DZUpperNonCornerGhostBox()
{
    return {{mhdGhostWidth, mhdGhostWidth, mhdGhostWidth + nCellsMHDZ3D},
            {mhdGhostWidth + nCellsMHDX3D - 1u, mhdGhostWidth + nCellsMHDY3D - 1u,
             2u * mhdGhostWidth + nCellsMHDZ3D - 1u}};
}


// ─── 1D FixedPressureOutflow ──────────────────────────────────────────────────

struct FixedPressureOutflowBC1D : testing::Test
{
    static constexpr double gamma   = 5.0 / 3.0;
    static constexpr double rho_val = 2.0;
    static constexpr double vx_val  = 1.0;
    static constexpr double vy_val  = 1.0;
    static constexpr double vz_val  = 1.0;
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

    GridLayoutMHD1D layout{{0.1}, {nCellsMHD}, {0.0}};

    GridMHD1D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD1D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD1D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<1> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<1> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<1> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<1>& rhoField{*(&rhoGrid)};
    FieldMHD<1>& PField{*(&PGrid)};
    FieldMHD<1>& EtotField{*(&EtotGrid)};

    FixedPressureOutflowBC1D()
    {
        auto fill_scalar = [&](auto& f, MHDQuantity::Scalar qty, double interior_val) {
            for (std::uint32_t i = 0; i < f.shape()[0]; ++i)
                f(i) = -999.0;
            std::uint32_t ps = layout.physicalStartIndex(qty, Direction::X);
            std::uint32_t pe = layout.physicalEndIndex(qty, Direction::X);
            for (std::uint32_t i = ps; i <= pe; ++i)
                f(i) = interior_val;
        };

        fill_scalar(rhoField, MHDQuantity::Scalar::rho, rho_val);
        fill_scalar(PField, MHDQuantity::Scalar::P, P_val);
        fill_scalar(EtotField, MHDQuantity::Scalar::Etot, etot_val);

        fill_scalar(rhoV[0], MHDQuantity::Scalar::rhoVx, rho_val * vx_val);
        fill_scalar(rhoV[1], MHDQuantity::Scalar::rhoVy, rho_val * vy_val);
        fill_scalar(rhoV[2], MHDQuantity::Scalar::rhoVz, rho_val * vz_val);

        fill_scalar(Bvec[0], MHDQuantity::Scalar::Bx, Bx_val);
        fill_scalar(Bvec[1], MHDQuantity::Scalar::By, By_val);
        fill_scalar(Bvec[2], MHDQuantity::Scalar::Bz, Bz_val);
    }
};

TEST_F(FixedPressureOutflowBC1D, DirichletPressureGhostEtotMatchesExpected)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
    auto P_bc
        = std::make_shared<FieldDirichletBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>(P_fixed);
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto B_bc = std::make_shared<
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhdUpperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto etotQty     = MHDQuantity::Scalar::Etot;
    std::uint32_t ps = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe = layout.physicalEndIndex(etotQty, Direction::X);

    for (std::uint32_t g = 0; g < mhdGhostWidth; ++g)
    {
        EXPECT_NEAR(EtotField(ps - 1u - g), etot_ghost, 1e-12) << "lower ghost g=" << g;
        EXPECT_NEAR(EtotField(pe + 1u + g), etot_ghost, 1e-12) << "upper ghost g=" << g;
    }
}

TEST_F(FixedPressureOutflowBC1D, InteriorEtotUnchangedAfterBC)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>();
    auto P_bc
        = std::make_shared<FieldDirichletBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D>>(P_fixed);
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto B_bc = std::make_shared<
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<1>, GridLayoutMHD1D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<1>, GridLayoutMHD1D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhdLowerGhostCellBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhdUpperGhostCellBox(), layout, makeCtx(acc, 0.0));

    auto etotQty     = MHDQuantity::Scalar::Etot;
    std::uint32_t ps = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe = layout.physicalEndIndex(etotQty, Direction::X);

    for (std::uint32_t i = ps; i <= pe; ++i)
        EXPECT_DOUBLE_EQ(EtotField(i), etot_val) << "interior index i=" << i;
}


// ─── 2D FixedPressureOutflow ──────────────────────────────────────────────────

struct FixedPressureOutflowBC2D : testing::Test
{
    static constexpr double gamma   = 5.0 / 3.0;
    static constexpr double rho_val = 2.0;
    static constexpr double vx_val  = 1.0;
    static constexpr double vy_val  = 1.0;
    static constexpr double vz_val  = 1.0;
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

    GridLayoutMHD2D layout{{0.1, 0.1}, {nCellsMHDX2D, nCellsMHDY2D}, {0.0, 0.0}};

    GridMHD2D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD2D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD2D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<2> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<2> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<2> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<2>& rhoField{*(&rhoGrid)};
    FieldMHD<2>& PField{*(&PGrid)};
    FieldMHD<2>& EtotField{*(&EtotGrid)};

    FixedPressureOutflowBC2D()
    {
        auto fill_all = [&](auto& f, double val) {
            auto [nx, ny] = f.shape();
            for (std::uint32_t ix = 0; ix < nx; ++ix)
                for (std::uint32_t iy = 0; iy < ny; ++iy)
                    f(ix, iy) = val;
        };

        fill_all(rhoField, rho_val);
        fill_all(PField, P_val);
        fill_all(EtotField, etot_val);

        fill_all(rhoV[0], rho_val * vx_val);
        fill_all(rhoV[1], rho_val * vy_val);
        fill_all(rhoV[2], rho_val * vz_val);

        fill_all(Bvec[0], Bx_val);
        fill_all(Bvec[1], By_val);
        fill_all(Bvec[2], Bz_val);
    }
};

TEST_F(FixedPressureOutflowBC2D, DirichletPressureGhostEtotMatchesExpected)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
    auto P_bc
        = std::make_shared<FieldDirichletBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>(P_fixed);
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto B_bc = std::make_shared<
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd2DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd2DYUpperGhostBox(), layout, makeCtx(acc, 0.0));

    for (auto const& idx : fp2DXLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "XLower ghost (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : fp2DXUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "XUpper ghost (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : fp2DYLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YLower ghost (" << idx[0] << "," << idx[1] << ")";
    for (auto const& idx : fp2DYUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YUpper ghost (" << idx[0] << "," << idx[1] << ")";
}

TEST_F(FixedPressureOutflowBC2D, InteriorEtotUnchangedAfterBC)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>();
    auto P_bc
        = std::make_shared<FieldDirichletBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D>>(P_fixed);
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto B_bc = std::make_shared<
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<2>, GridLayoutMHD2D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<2>, GridLayoutMHD2D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd2DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd2DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd2DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd2DYUpperGhostBox(), layout, makeCtx(acc, 0.0));

    auto etotQty       = MHDQuantity::Scalar::Etot;
    std::uint32_t ps_x = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe_x = layout.physicalEndIndex(etotQty, Direction::X);
    std::uint32_t ps_y = layout.physicalStartIndex(etotQty, Direction::Y);
    std::uint32_t pe_y = layout.physicalEndIndex(etotQty, Direction::Y);

    for (std::uint32_t ix = ps_x; ix <= pe_x; ++ix)
        for (std::uint32_t iy = ps_y; iy <= pe_y; ++iy)
            EXPECT_DOUBLE_EQ(EtotField(ix, iy), etot_val)
                << "interior index (" << ix << "," << iy << ")";
}


// ─── 3D FixedPressureOutflow ──────────────────────────────────────────────────

struct FixedPressureOutflowBC3D : testing::Test
{
    static constexpr double gamma   = 5.0 / 3.0;
    static constexpr double rho_val = 2.0;
    static constexpr double vx_val  = 1.0;
    static constexpr double vy_val  = 1.0;
    static constexpr double vz_val  = 1.0;
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

    GridLayoutMHD3D layout{
        {0.1, 0.1, 0.1}, {nCellsMHDX3D, nCellsMHDY3D, nCellsMHDZ3D}, {0.0, 0.0, 0.0}};

    GridMHD3D rhoGrid{"rho", MHDQuantity::Scalar::rho, layout.allocSize(MHDQuantity::Scalar::rho)};
    GridMHD3D PGrid{"P", MHDQuantity::Scalar::P, layout.allocSize(MHDQuantity::Scalar::P)};
    GridMHD3D EtotGrid{"Etot", MHDQuantity::Scalar::Etot,
                       layout.allocSize(MHDQuantity::Scalar::Etot)};

    UsableVecFieldMHD<3> rhoV{"rhoV", layout, MHDQuantity::Vector::rhoV};
    UsableVecFieldMHD<3> Bvec{"B", layout, MHDQuantity::Vector::B};

    MHDPatchFieldAccessorTest<3> acc{rhoGrid, PGrid, EtotGrid, rhoV, Bvec};

    FieldMHD<3>& rhoField{*(&rhoGrid)};
    FieldMHD<3>& PField{*(&PGrid)};
    FieldMHD<3>& EtotField{*(&EtotGrid)};

    FixedPressureOutflowBC3D()
    {
        auto fill_all = [&](auto& f, double val) {
            auto [nx, ny, nz] = f.shape();
            for (std::uint32_t ix = 0; ix < nx; ++ix)
                for (std::uint32_t iy = 0; iy < ny; ++iy)
                    for (std::uint32_t iz = 0; iz < nz; ++iz)
                        f(ix, iy, iz) = val;
        };

        fill_all(rhoField, rho_val);
        fill_all(PField, P_val);
        fill_all(EtotField, etot_val);

        fill_all(rhoV[0], rho_val * vx_val);
        fill_all(rhoV[1], rho_val * vy_val);
        fill_all(rhoV[2], rho_val * vz_val);

        fill_all(Bvec[0], Bx_val);
        fill_all(Bvec[1], By_val);
        fill_all(Bvec[2], Bz_val);
    }
};

TEST_F(FixedPressureOutflowBC3D, DirichletPressureGhostEtotMatchesExpected)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
    auto P_bc
        = std::make_shared<FieldDirichletBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>(P_fixed);
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto B_bc = std::make_shared<
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd3DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd3DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd3DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd3DYUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZLower, mhd3DZLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZUpper, mhd3DZUpperGhostBox(), layout, makeCtx(acc, 0.0));

    for (auto const& idx : fp3DXLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "XLower ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : fp3DXUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "XUpper ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : fp3DYLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YLower ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : fp3DYUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "YUpper ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : fp3DZLowerNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "ZLower ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
    for (auto const& idx : fp3DZUpperNonCornerGhostBox())
        EXPECT_NEAR(EtotField(idx), etot_ghost, 1e-12)
            << "ZUpper ghost (" << idx[0] << "," << idx[1] << "," << idx[2] << ")";
}

TEST_F(FixedPressureOutflowBC3D, InteriorEtotUnchangedAfterBC)
{
    auto rho_bc = std::make_shared<FieldNeumannBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>();
    auto P_bc
        = std::make_shared<FieldDirichletBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D>>(P_fixed);
    auto rhoV_bc
        = std::make_shared<FieldNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto B_bc = std::make_shared<
        FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldMHD<3>, GridLayoutMHD3D>>();
    auto thermo = std::make_shared<IdealGasThermo>(gamma);

    FieldTotalEnergyFromPressureBoundaryCondition<FieldMHD<3>, GridLayoutMHD3D> bc{
        rho_bc, rhoV_bc, B_bc, P_bc, thermo};

    bc.apply(EtotField, BoundaryLocation::XLower, mhd3DXLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::XUpper, mhd3DXUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YLower, mhd3DYLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::YUpper, mhd3DYUpperGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZLower, mhd3DZLowerGhostBox(), layout, makeCtx(acc, 0.0));
    bc.apply(EtotField, BoundaryLocation::ZUpper, mhd3DZUpperGhostBox(), layout, makeCtx(acc, 0.0));

    auto etotQty       = MHDQuantity::Scalar::Etot;
    std::uint32_t ps_x = layout.physicalStartIndex(etotQty, Direction::X);
    std::uint32_t pe_x = layout.physicalEndIndex(etotQty, Direction::X);
    std::uint32_t ps_y = layout.physicalStartIndex(etotQty, Direction::Y);
    std::uint32_t pe_y = layout.physicalEndIndex(etotQty, Direction::Y);
    std::uint32_t ps_z = layout.physicalStartIndex(etotQty, Direction::Z);
    std::uint32_t pe_z = layout.physicalEndIndex(etotQty, Direction::Z);

    for (std::uint32_t ix = ps_x; ix <= pe_x; ++ix)
        for (std::uint32_t iy = ps_y; iy <= pe_y; ++iy)
            for (std::uint32_t iz = ps_z; iz <= pe_z; ++iz)
                EXPECT_DOUBLE_EQ(EtotField(ix, iy, iz), etot_val)
                    << "interior index (" << ix << "," << iy << "," << iz << ")";
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
