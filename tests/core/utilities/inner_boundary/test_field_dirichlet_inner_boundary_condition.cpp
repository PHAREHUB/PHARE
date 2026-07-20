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
#include "core/inner_boundary/field_dirichlet_inner_boundary_condition.hpp"
#include "core/utilities/box/box.hpp"

namespace
{
constexpr double eps = 1e-12;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Classifier     = PHARE::core::InnerBoundaryMeshClassifier<2, GridLayout, PHARE::core::MHDQuantity>;
using MeshData       = PHARE::core::InnerBoundaryMeshData<2, PHARE::core::MHDQuantity>;
using ScalarField    = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;

struct DummyState
{
};

constexpr std::array<PHARE::core::QtyCentering, 2> kCellC
    = {PHARE::core::QtyCentering::dual, PHARE::core::QtyCentering::dual};

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

struct DirichletBCFixture
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    MeshDataBuffers buffers{layout};

    DirichletBCFixture()
    {
        Classifier::Overrides ov;
        ov.cut_eps      = 1e-12;
        ov.inactive_eps = 0.0;
        auto classifier = Classifier::withDefaults(plane, layout, ov);
        classifier(layout, buffers.tags);
    }
};

} // namespace

TEST(FieldDirichletInnerBoundaryCondition, constantFieldMatchingBoundaryValueIsUnchanged)
{
    DirichletBCFixture fix;
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

    PHARE::core::FieldDirichletInnerBoundaryCondition<ScalarField, GridLayout, DummyState> bc{C};
    DummyState state;
    bc.apply(field, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            EXPECT_NEAR(field(i, j), C, eps) << "cell (" << i << ", " << j << ") changed";
}

TEST(FieldDirichletInnerBoundaryCondition, ghostCellReceivesExtrapolatedBoundaryValue)
{
    DirichletBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    constexpr double boundaryValue = 5.0;

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

    PHARE::core::FieldDirichletInnerBoundaryCondition<ScalarField, GridLayout, DummyState> bc{
        boundaryValue};
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
                << ") must not be touched by the BC";
            continue;
        }

        foundInPatch          = true;
        double const expected = 2.0 * boundaryValue - (g.mirrorPoint[0] + g.mirrorPoint[1]);
        EXPECT_NEAR(field(g.index), expected, 1e-10)
            << "ghost at (" << g.index[0] << ", " << g.index[1] << ")";
    }
    EXPECT_TRUE(foundInPatch) << "at least one in-patch ghost must exist";
}

TEST(FieldDirichletInnerBoundaryCondition, fillsNonInterpolableGhostsReflectsExtrapolationMode)
{
    using BC = PHARE::core::FieldDirichletInnerBoundaryCondition<ScalarField, GridLayout, DummyState>;
    BC constantBC{0.0, BC::ExtrapolationType::Constant};
    BC linearBC{0.0, BC::ExtrapolationType::Linear};
    EXPECT_TRUE(constantBC.fillsNonInterpolableGhosts());
    EXPECT_FALSE(linearBC.fillsNonInterpolableGhosts());
}

TEST(FieldDirichletInnerBoundaryCondition, constantModeFillsEveryGhostIncludingNonInterpolable)
{
    DirichletBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    constexpr double boundaryValue = 3.0;
    constexpr double sentinel      = -999.0;

    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)};
    ScalarField field{"rho", PHARE::core::MHDQuantity::Scalar::CellCentered,
                      storage.data(), storage.shape()};

    for (auto const& g : meshData.getGhostDataFromCentering(kCellC))
        field(g.index) = sentinel;

    using BC = PHARE::core::FieldDirichletInnerBoundaryCondition<ScalarField, GridLayout, DummyState>;
    BC bc{boundaryValue, BC::ExtrapolationType::Constant};
    DummyState state;
    bc.apply(field, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    auto const& ghostCells = meshData.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghostCells.empty());

    bool sawNonInterpolable = false;
    for (auto const& g : ghostCells)
    {
        if (!g.mirrorIsInterpolable)
            sawNonInterpolable = true;
        // Constant mode must fill *every* ghost, interpolable or not.
        EXPECT_NEAR(field(g.index), boundaryValue, eps)
            << "ghost at (" << g.index[0] << ", " << g.index[1]
            << ") mirrorIsInterpolable=" << g.mirrorIsInterpolable << " not filled in constant mode";
    }
    EXPECT_TRUE(sawNonInterpolable)
        << "fixture should expose at least one non-interpolable ghost to exercise the fix";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
