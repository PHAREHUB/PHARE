#include "gtest/gtest.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/inner_boundary/field_neumann_inner_boundary_condition.hpp"
#include "core/utilities/box/box.hpp"

namespace
{
constexpr double eps = 1e-12;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Classifier     = PHARE::core::InnerBoundaryMeshClassifier<2, GridLayout, PHARE::core::MHDQuantity>;
using MeshData       = PHARE::core::InnerBoundaryMeshData<2, PHARE::core::MHDQuantity>;
using ScalarField    = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;

/// Minimal physical-state stub; the Neumann BC does not use the state argument.
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
 * @brief Standard 2D plane-boundary test fixture.
 *
 * Plane at x = 0, outward normal (1, 0). Grid: 4 × 2 physical cells with
 * cell width dx = dy = 1. AMR box {{-2,0},{1,1}} → physical cell centres at
 * x = {-1.5, -0.5, 0.5, 1.5}. The ghost cell at physical[0,0] has its mirror
 * at (1.5, 0.5) which lies inside the fluid region.
 */
struct NeumannBCFixture
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    MeshDataBuffers buffers{layout};

    NeumannBCFixture()
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
// Tests
// ---------------------------------------------------------------------------

/**
 * @brief A spatially constant field must be unchanged after applying a Neumann BC.
 *
 * For a homogeneous Neumann condition (∂f/∂n = 0), the BC sets
 * f(ghost) = f(mirrorPoint).  When every mesh node already carries the same
 * constant value @c C, the interpolated value at any mirror point is also @c C,
 * so the ghost values remain @c C.
 */
TEST(FieldNeumannInnerBoundaryCondition, constantFieldIsUnchanged)
{
    NeumannBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    constexpr double C = 7.0;

    // Allocate a cell-centred scalar field and fill every node (including AMR
    // ghost cells) with the constant value C.
    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)};
    ScalarField field{"rho", PHARE::core::MHDQuantity::Scalar::CellCentered,
                      storage.data(), storage.shape()};

    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            field(i, j) = C;

    PHARE::core::FieldNeumannInnerBoundaryCondition<ScalarField, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(field, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    // Every cell must still hold the constant value.
    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
            EXPECT_NEAR(field(i, j), C, eps) << "cell (" << i << ", " << j << ") changed";
}

/**
 * @brief The Neumann BC correctly copies the mirror-point value into the ghost cell.
 *
 * We set the fluid region to a non-trivial linear profile f(x,y) = x + y and
 * the ghost region to zero.  After applying the Neumann BC the ghost at
 * physical[0,0] (centre (-1.5, 0.5)) must receive the value interpolated at
 * its mirror point (1.5, 0.5), i.e. f = 1.5 + 0.5 = 2.0.
 */
TEST(FieldNeumannInnerBoundaryCondition, ghostCellReceivesMirrorPointValue)
{
    NeumannBCFixture fix;
    auto const& layout   = fix.layout;
    auto const& meshData = fix.buffers.tags;

    // Allocate a cell-centred scalar field.
    PHARE::core::NdArrayVector<2, double> storage{
        layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)};
    ScalarField field{"rho", PHARE::core::MHDQuantity::Scalar::CellCentered,
                      storage.data(), storage.shape()};

    // Fill every node with f(x,y) = x + y using fieldNodeCoordinates.
    for (auto i = 0u; i < field.shape()[0]; ++i)
        for (auto j = 0u; j < field.shape()[1]; ++j)
        {
            auto amr_pos = layout.localToAMR(PHARE::core::Point<std::uint32_t, 2>{i, j});
            auto amr_ij  = PHARE::core::Point<int, 2>{static_cast<int>(amr_pos[0]),
                                                      static_cast<int>(amr_pos[1])};
            auto pos     = layout.fieldNodeCoordinates(field, amr_ij);
            field(i, j)  = pos[0] + pos[1];
        }

    // Zero the ghost cells so the change is observable.
    for (auto const& g : meshData.getGhostDataFromCentering(kCellC))
        field(g.index) = 0.0;

    PHARE::core::FieldNeumannInnerBoundaryCondition<ScalarField, GridLayout, DummyState> bc;
    DummyState state;
    bc.apply(field, layout, meshData, PHARE::core::InnerBCContext<DummyState>{state, state, 0.0});

    // Only in-patch ghosts should be updated; out-of-patch ones stay at 0 (filled later by AMR).
    auto const& ghost_cells = meshData.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghost_cells.empty());

    bool found_in_patch = false;
    for (auto const& g : ghost_cells)
    {
        if (!g.mirrorIsInterpolable)
        {
            EXPECT_NEAR(field(g.index), 0.0, eps)
                << "out-of-patch ghost at (" << g.index[0] << ", " << g.index[1]
                << ") must not be touched by the BC";
            continue;
        }
        found_in_patch          = true;
        double const expected   = g.mirrorPoint[0] + g.mirrorPoint[1];
        EXPECT_NEAR(field(g.index), expected, 1e-10)
            << "ghost at (" << g.index[0] << ", " << g.index[1] << ")";
    }
    EXPECT_TRUE(found_in_patch) << "at least one in-patch ghost must exist";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
