#include "gtest/gtest.h"

#include <cstdint>
#include <vector>

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/inner_boundary/field_symmetric_inner_boundary_condition.hpp"
#include "core/numerics/interpolator/field_at_point.hpp"
#include "core/utilities/box/box.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures_mhd.hpp"

namespace
{
constexpr double eps = 1e-10;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Classifier     = PHARE::core::InnerBoundaryMeshClassifier<2, GridLayout, PHARE::core::MHDQuantity>;
using MeshData       = PHARE::core::InnerBoundaryMeshData<2, PHARE::core::MHDQuantity>;
using ScalarField    = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using VecFieldMHD2   = PHARE::core::UsableVecFieldMHD<2>;

struct DummyState
{
};

constexpr std::array<PHARE::core::QtyCentering, 2> kCellC
    = {PHARE::core::QtyCentering::dual, PHARE::core::QtyCentering::dual};

/**
 * @brief Allocates and wires all mesh-data buffers needed by InnerBoundaryMeshData.
 */
struct MeshDataBuffers
{
    static constexpr char const* BOUNDARY_NAME = "test";

    explicit MeshDataBuffers(GridLayout const& layout)
        : phi_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::NodeCentered)}
        , tags{BOUNDARY_NAME}
    {
        ScalarField phi_field{std::string(BOUNDARY_NAME) + "_signed_distance",
                              PHARE::core::MHDQuantity::Scalar::NodeCentered,
                              phi_storage.data(), phi_storage.shape()};
        tags.signedDistanceAtNodes.setBuffer(&phi_field);

        elem_storages.reserve(MeshData::num_elem_types);
        for (std::size_t i = 0; i < MeshData::num_elem_types; ++i)
        {
            auto const c   = MeshData::idxToCentering(i);
            auto const qty = MeshData::scalarFromCentering(c);
            elem_storages.emplace_back(layout.allocSize(qty));
            ScalarField tmp{tags.elemStatus[i].name(), qty,
                            elem_storages[i].data(), elem_storages[i].shape()};
            tags.elemStatus[i].setBuffer(&tmp);
        }
        tags.ghostElemsData._data    = &ghost_array;
        tags.degradedElemsData._data = &degraded_array;
    }

    PHARE::core::NdArrayVector<2, double> phi_storage;
    std::vector<PHARE::core::NdArrayVector<2, double>> elem_storages;
    PHARE::core::GhostElemPack<2>::ghost_elem_array_type ghost_array{};
    PHARE::core::GhostElemPack<2>::ghost_elem_array_type degraded_array{};
    MeshData tags;
};

/**
 * @brief Plane at x = 0, outward normal (1, 0). Grid: 4 × 2 physical cells,
 *        dx = dy = 1. AMR box {{-2,0},{1,1}} → cell centres at
 *        x = {-1.5, -0.5, 0.5, 1.5}. Ghost cells on the left (x < 0) have
 *        their mirrors on the right (x > 0); specifically the in-patch ghost
 *        at local[0,j] has its mirror at (1.5, j_center).
 */
struct SymmetricBCFixture
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    MeshDataBuffers buffers{layout};

    SymmetricBCFixture()
    {
        Classifier::Overrides ov;
        ov.cut_eps      = 1e-12;
        ov.inactive_eps = 0.0;
        auto classifier = Classifier::withDefaults(plane, layout, ov);
        classifier(layout, buffers.tags);
    }
};

} // namespace

// ---------------------------------------------------------------------------
// Scalar tests (falls back to Neumann)
// ---------------------------------------------------------------------------

/**
 * @brief A spatially constant scalar field is unchanged after applying the symmetric BC.
 *
 * The scalar path falls back to Neumann (∂f/∂n = 0), which copies the mirror
 * value into the ghost. When the field is constant the mirror value equals the
 * ghost value, so nothing changes.
 */
TEST(FieldSymmetricInnerBoundaryCondition, scalarConstantFieldIsUnchanged)
{
    SymmetricBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    constexpr double C = 7.0;

    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)};
    ScalarField field{"rho", PHARE::core::MHDQuantity::Scalar::CellCentered,
                      storage.data(), storage.shape()};

    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            field(i, j) = C;

    PHARE::core::FieldSymmetricInnerBoundaryCondition<ScalarField, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(field, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            EXPECT_NEAR(field(i, j), C, eps) << "cell (" << i << ", " << j << ") changed";
}

/**
 * @brief Ghost scalar cell receives the mirror-point value (Neumann fallback).
 */
TEST(FieldSymmetricInnerBoundaryCondition, scalarGhostCellReceivesMirrorPointValue)
{
    SymmetricBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)};
    ScalarField field{"rho", PHARE::core::MHDQuantity::Scalar::CellCentered,
                      storage.data(), storage.shape()};

    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
        {
            auto amr_pos = layout.localToAMR(PHARE::core::Point<std::uint32_t, 2>{i, j});
            auto amr_ij  = PHARE::core::Point<int, 2>{static_cast<int>(amr_pos[0]),
                                                      static_cast<int>(amr_pos[1])};
            auto pos     = layout.fieldNodeCoordinates(field, amr_ij);
            field(i, j)  = pos[0] + pos[1];
        }

    for (auto const& g : meshData.getGhostDataFromCentering(kCellC))
        field(g.index) = 0.0;

    PHARE::core::FieldSymmetricInnerBoundaryCondition<ScalarField, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(field, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    auto const& ghostCells = meshData.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghostCells.empty());

    bool foundInPatch = false;
    for (auto const& g : ghostCells)
    {
        if (!g.mirrorIsInterpolable)
        {
            EXPECT_NEAR(field(g.index), 0.0, eps)
                << "out-of-patch ghost at (" << g.index[0] << ", " << g.index[1]
                << ") must not be touched";
            continue;
        }
        foundInPatch         = true;
        double const expected = g.mirrorPoint[0] + g.mirrorPoint[1];
        EXPECT_NEAR(field(g.index), expected, eps)
            << "ghost at (" << g.index[0] << ", " << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

// ---------------------------------------------------------------------------
// Vector tests
// ---------------------------------------------------------------------------

/**
 * @brief For a constant vector field (Cx, Cy, Cz), the symmetric BC on a plane
 *        with normal n=(1,0) reverses the x-component and preserves y and z.
 *
 * All 3 components are CellCentered so they share the same ghost-cell list.
 *
 * For V = (Cx, Cy, Cz) and normal n=(1,0,0):
 *   V_n = (Cx, 0, 0),  V_t = (0, Cy, Cz)
 *   V_reflected = V_t - V_n = (-Cx, Cy, Cz)
 */
TEST(FieldSymmetricInnerBoundaryCondition, vectorConstantField_normalComponentFlips)
{
    SymmetricBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    constexpr double Cx = 3.0, Cy = 4.0, Cz = 5.0;

    VecFieldMHD2 V{"V", layout, PHARE::core::MHDQuantity::Vector::CellCentered};
    auto& Vx = V.getComponent(PHARE::core::Component::X);
    auto& Vy = V.getComponent(PHARE::core::Component::Y);
    auto& Vz = V.getComponent(PHARE::core::Component::Z);

    for (auto i = 0u; i < Vx.shape()[0]; ++i)
        for (auto j = 0u; j < Vx.shape()[1]; ++j)
        {
            Vx(i, j) = Cx;
            Vy(i, j) = Cy;
            Vz(i, j) = Cz;
        }

    PHARE::core::FieldSymmetricInnerBoundaryCondition<VecFieldMHD2, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(V, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    auto const& ghostCells = meshData.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghostCells.empty());

    bool foundInPatch = false;
    for (auto const& g : ghostCells)
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch = true;

        EXPECT_NEAR(Vx(g.index), -Cx, eps)
            << "Vx ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(Vy(g.index), Cy, eps)
            << "Vy ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(Vz(g.index), Cz, eps)
            << "Vz ghost at (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

/**
 * @brief For a vector with no normal component, the symmetric BC leaves the
 *        field unchanged.
 *
 * If V = (0, Vy, Vz) and n=(1,0,0), the normal component is zero, so the
 * reflection V_t - V_n = V_t = V: the ghost value equals the mirror value.
 */
TEST(FieldSymmetricInnerBoundaryCondition, vectorPurelyTangentialField_isUnchanged)
{
    SymmetricBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    constexpr double Cy = 2.0, Cz = -3.0;

    VecFieldMHD2 V{"V", layout, PHARE::core::MHDQuantity::Vector::CellCentered};
    auto& Vx = V.getComponent(PHARE::core::Component::X);
    auto& Vy = V.getComponent(PHARE::core::Component::Y);
    auto& Vz = V.getComponent(PHARE::core::Component::Z);

    for (auto i = 0u; i < Vx.shape()[0]; ++i)
        for (auto j = 0u; j < Vx.shape()[1]; ++j)
        {
            Vx(i, j) = 0.0;
            Vy(i, j) = Cy;
            Vz(i, j) = Cz;
        }

    PHARE::core::FieldSymmetricInnerBoundaryCondition<VecFieldMHD2, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(V, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    auto const& ghostCells = meshData.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghostCells.empty());

    bool foundInPatch = false;
    for (auto const& g : ghostCells)
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch = true;

        EXPECT_NEAR(Vx(g.index), 0.0, eps)
            << "Vx ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(Vy(g.index), Cy, eps)
            << "Vy ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(Vz(g.index), Cz, eps)
            << "Vz ghost at (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

/**
 * @brief For a linear vector field, the ghost receives the correctly reflected
 *        interpolated mirror value.
 *
 * Field: Vx(x,y) = x+y, Vy(x,y) = x-y, Vz = 1.
 * Boundary normal n = (1, 0) → extended to (1, 0, 0) in 3D.
 *
 * For an in-patch ghost whose mirror is at (1.5, y_m):
 *   V_mirror = (1.5 + y_m,  1.5 - y_m,  1)
 *   V_n      = (1.5 + y_m,  0,           0)
 *   V_t      = (0,           1.5 - y_m,  1)
 *   V_refl   = (-( 1.5+y_m),  1.5-y_m,   1)
 */
TEST(FieldSymmetricInnerBoundaryCondition, vectorLinearField_correctReflection)
{
    SymmetricBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    VecFieldMHD2 V{"V", layout, PHARE::core::MHDQuantity::Vector::CellCentered};
    auto& Vx = V.getComponent(PHARE::core::Component::X);
    auto& Vy = V.getComponent(PHARE::core::Component::Y);
    auto& Vz = V.getComponent(PHARE::core::Component::Z);

    for (auto i = 0u; i < Vx.shape()[0]; ++i)
        for (auto j = 0u; j < Vx.shape()[1]; ++j)
        {
            auto amr_pos = layout.localToAMR(PHARE::core::Point<std::uint32_t, 2>{i, j});
            auto amr_ij  = PHARE::core::Point<int, 2>{static_cast<int>(amr_pos[0]),
                                                      static_cast<int>(amr_pos[1])};
            auto pos     = layout.fieldNodeCoordinates(Vx, amr_ij);
            Vx(i, j) = pos[0] + pos[1];
            Vy(i, j) = pos[0] - pos[1];
            Vz(i, j) = 1.0;
        }

    // Zero ghost cells so changes are observable.
    for (auto const& g : meshData.getGhostDataFromCentering(kCellC))
    {
        Vx(g.index) = 0.0;
        Vy(g.index) = 0.0;
        Vz(g.index) = 0.0;
    }

    PHARE::core::FieldSymmetricInnerBoundaryCondition<VecFieldMHD2, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(V, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    auto const& ghostCells = meshData.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghostCells.empty());

    bool foundInPatch = false;
    for (auto const& g : ghostCells)
    {
        if (!g.mirrorIsInterpolable)
            continue;
        foundInPatch = true;

        double const xm = g.mirrorPoint[0];
        double const ym = g.mirrorPoint[1];
        // V_mirror = (xm+ym, xm-ym, 1); n=(1,0,0)
        // V_refl   = (-(xm+ym), xm-ym, 1)
        double const exp_vx = -(xm + ym);
        double const exp_vy =  (xm - ym);
        double const exp_vz =  1.0;

        EXPECT_NEAR(Vx(g.index), exp_vx, eps)
            << "Vx ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(Vy(g.index), exp_vy, eps)
            << "Vy ghost at (" << g.index[0] << "," << g.index[1] << ")";
        EXPECT_NEAR(Vz(g.index), exp_vz, eps)
            << "Vz ghost at (" << g.index[0] << "," << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

/**
 * @brief Regression: when a ghost's mirror is interpolable for component i but NOT for some
 *        sibling component j (different centering), the BC must skip that ghost rather than
 *        read OOB while interpolating component j.
 *
 * See the matching test in the antisymmetric BC file for the rationale.
 */
TEST(FieldSymmetricInnerBoundaryCondition, vectorMixedCenterings_skipGhostWhenSiblingNotInterpolable)
{
    using PHARE::core::Component;
    using PHARE::core::FieldAtPoint;
    using PHARE::core::MHDQuantity;
    using PHARE::core::Point;
    using PHARE::core::QtyCentering;

    SymmetricBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    VecFieldMHD2 B{"B", layout, MHDQuantity::Vector::B};
    auto& Bx = B.getComponent(Component::X);
    auto& By = B.getComponent(Component::Y);
    auto& Bz = B.getComponent(Component::Z);

    std::array<std::array<QtyCentering, 2>, 3> const compCentering{
        {GridLayout::centering(MHDQuantity::Scalar::Bx),
         GridLayout::centering(MHDQuantity::Scalar::By),
         GridLayout::centering(MHDQuantity::Scalar::Bz)}};

    struct Trigger
    {
        std::size_t i;
        Point<std::uint32_t, 2> idx;
    };
    std::vector<Trigger> triggers;
    for (std::size_t i = 0; i < 3; ++i)
    {
        for (auto const& g : meshData.getGhostDataFromCentering(compCentering[i]))
        {
            if (!g.mirrorIsInterpolable)
                continue;
            for (std::size_t j = 0; j < 3; ++j)
            {
                if (j == i)
                    continue;
                if (!FieldAtPoint<2, 1>::pointIsInterpolable(layout, g.mirrorPoint,
                                                             compCentering[j]))
                {
                    triggers.push_back({i, g.index});
                    break;
                }
            }
        }
    }

    ASSERT_FALSE(triggers.empty())
        << "fixture geometry produces no ghost where one component is interpolable but a "
           "sibling is not; cannot exercise the sibling-check regression";

    constexpr double sentinel = 12345.6789;
    for (auto i = 0u; i < Bx.shape()[0]; ++i)
        for (auto j = 0u; j < Bx.shape()[1]; ++j)
            Bx(i, j) = 1.0;
    for (auto i = 0u; i < By.shape()[0]; ++i)
        for (auto j = 0u; j < By.shape()[1]; ++j)
            By(i, j) = 2.0;
    for (auto i = 0u; i < Bz.shape()[0]; ++i)
        for (auto j = 0u; j < Bz.shape()[1]; ++j)
            Bz(i, j) = 3.0;

    for (auto const& t : triggers)
    {
        if (t.i == 0)
            Bx(t.idx) = sentinel;
        else if (t.i == 1)
            By(t.idx) = sentinel;
        else
            Bz(t.idx) = sentinel;
    }

    PHARE::core::FieldSymmetricInnerBoundaryCondition<VecFieldMHD2, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(B, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    for (auto const& t : triggers)
    {
        double const got = (t.i == 0) ? Bx(t.idx) : (t.i == 1) ? By(t.idx) : Bz(t.idx);
        EXPECT_NEAR(got, sentinel, eps)
            << "component i=" << t.i << " at local (" << t.idx[0] << "," << t.idx[1]
            << ") was overwritten despite a sibling being non-interpolable";
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
