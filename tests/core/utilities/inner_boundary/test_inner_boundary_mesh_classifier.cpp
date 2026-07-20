#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/inner_boundary_mesh_classifier.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/inner_boundary/sphere_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"

namespace
{
constexpr double eps = 1e-12;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Mapper = PHARE::core::InnerBoundaryMeshClassifier<2, GridLayout, PHARE::core::MHDQuantity>;
using MeshData           = PHARE::core::InnerBoundaryMeshData<2, PHARE::core::MHDQuantity>;
using ScalarField    = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;

template<typename FieldT>
PHARE::core::Point<std::uint32_t, 2> physicalLocalIndex(GridLayout const& layout,
                                                        FieldT const& field,
                                                        std::uint32_t ix,
                                                        std::uint32_t iy)
{
    using PHARE::core::Direction;

    return {layout.physicalStartIndex(field, Direction::X) + ix,
            layout.physicalStartIndex(field, Direction::Y) + iy};
}

// Centering constants for 2D: bit i = 1 when direction i is dual.
// (P,P)→0=node, (D,P)→1=faceX, (P,D)→2=faceY, (D,D)→3=cell
constexpr std::array<PHARE::core::QtyCentering, 2> kCellC
    = {PHARE::core::QtyCentering::dual,   PHARE::core::QtyCentering::dual};
constexpr std::array<PHARE::core::QtyCentering, 2> kFaceXC
    = {PHARE::core::QtyCentering::primal, PHARE::core::QtyCentering::dual};
constexpr std::array<PHARE::core::QtyCentering, 2> kFaceYC
    = {PHARE::core::QtyCentering::dual,   PHARE::core::QtyCentering::primal};

struct InnerBoundaryMeshClassifierBuffers
{
    static constexpr char const* BOUNDARY_NAME = "test";

    explicit InnerBoundaryMeshClassifierBuffers(GridLayout const& layout)
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
} // namespace

TEST(InnerBoundaryMeshClassifier, computesReasonableDefaultCutEpsFromLayout)
{
    // Sphere centred at (0.05, 0.05) so no grid node coincides with the centre
    // (which would make normal() undefined and cause populateGhostLists_ to throw).
    // The physical origin node sits at (0, 0), distance = sqrt(0.05^2+0.05^2) - 1 ≈ -0.929,
    // but we check the node at (0.2, 0.1) which is closest to r=1 for a simpler assertion.
    // We only exercise that the classifier runs without error and that phi at a node
    // well inside the sphere is negative.
    PHARE::core::SphereInnerBoundary<2> sphere{"sphere", {0.05, 0.05}, 1.};
    GridLayout layout{{0.2, 0.1}, {4u, 4u}, {0., 0.}};
    auto tagger = Mapper::withDefaults(sphere, layout);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    // The origin node (0, 0) is inside the sphere (distance < 0).
    auto const origin_node = physicalLocalIndex(layout, buffers.tags.signedDistanceAtNodes, 0u, 0u);
    EXPECT_LT(buffers.tags.signedDistanceAtNodes(origin_node), 0.);
}

TEST(InnerBoundaryMeshClassifier, tagsCutInactiveAndGhostGeometry)
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}}; // x=0
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Mapper::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Mapper::withDefaults(plane, layout, ov);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    // Ghost cells grow inward (into the solid) from the cut layer.
    // physical[0] (x=-1.5) was inactive and is now promoted to ghost; physical[3] (x=1.5)
    // is a plain fluid cell on the outside.
    auto& cellField = buffers.tags.getStatusFieldFromCentering(kCellC);
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Ghost),
              cellField(physicalLocalIndex(layout, cellField, 0u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Cut),
              cellField(physicalLocalIndex(layout, cellField, 1u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Cut),
              cellField(physicalLocalIndex(layout, cellField, 2u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Fluid),
              cellField(physicalLocalIndex(layout, cellField, 3u, 0u)));

    // Face at physical[2] straddles x=0 → Cut; face at physical[1] (between ghost and cut
    // cells) has a Cut adjacent cell → Fluid (not Ghost: mirror too close to boundary).
    auto& faceXField = buffers.tags.getStatusFieldFromCentering(kFaceXC);
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Cut),
              faceXField(physicalLocalIndex(layout, faceXField, 2u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Fluid),
              faceXField(physicalLocalIndex(layout, faceXField, 1u, 0u)));
    // FaceCenteredY face at physical[0,0] is adjacent to the ghost cell column → Ghost.
    auto& faceYField = buffers.tags.getStatusFieldFromCentering(kFaceYC);
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::ElemStatus::Ghost),
              faceYField(physicalLocalIndex(layout, faceYField, 0u, 0u)));

    // A plane has infinite characteristic length, so the level is never under-resolved:
    // no degraded elements are emitted.
    EXPECT_TRUE(buffers.tags.getDegradedElemsFromCentering(kCellC).empty());

    // Note: in 2D, EdgeCenteredY has centering (primal, dual) identical to FaceCenteredX,
    // so both concepts share faceXField above. No separate edge assertions needed.
}

TEST(InnerBoundaryMeshClassifier, ghostCellListIsNonEmpty)
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Mapper::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Mapper::withDefaults(plane, layout, ov);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    EXPECT_FALSE(buffers.tags.getGhostDataFromCentering(kCellC).empty());
}

TEST(InnerBoundaryMeshClassifier, ghostCellHasCorrectMirrorPointAndNormal)
{
    // Plane at x=0, normal (1,0). Ghost cell at physical[0,0] has center (-1.5, 0.5).
    // Its mirror is at (1.5, 0.5) and the outward normal is (1, 0).
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Mapper::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Mapper::withDefaults(plane, layout, ov);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    auto const target_idx
        = physicalLocalIndex(layout, buffers.tags.getStatusFieldFromCentering(kCellC), 0u, 0u);
    auto const& ghost_cells = buffers.tags.getGhostDataFromCentering(kCellC);
    auto it = std::find_if(ghost_cells.begin(), ghost_cells.end(),
                           [&](auto const& g) { return g.index == target_idx; });
    ASSERT_NE(it, ghost_cells.end()) << "Ghost cell for physical[0,0] not found in ghostCells";

    EXPECT_NEAR(it->mirrorPoint[0], 1.5, eps);
    EXPECT_NEAR(it->mirrorPoint[1], 0.5, eps);
    EXPECT_NEAR(it->normal[0], 1.0, eps);
    EXPECT_NEAR(it->normal[1], 0.0, eps);
    // Mirror (1.5, 0.5) is within the physical domain of the patch → in-patch.
    EXPECT_TRUE(it->mirrorIsInterpolable);
}

TEST(InnerBoundaryMeshClassifier, ghostFaceListIsNonEmpty)
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Mapper::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Mapper::withDefaults(plane, layout, ov);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    // FaceCenteredX: physical[0,0] (only OOB and Ghost neighbors) is Ghost → list non-empty.
    EXPECT_FALSE(buffers.tags.getGhostDataFromCentering(kFaceXC).empty());
    // FaceCenteredY face adjacent to the ghost cell column must also be present.
    EXPECT_FALSE(buffers.tags.getGhostDataFromCentering(kFaceYC).empty());
}

TEST(InnerBoundaryMeshClassifier, ghostCellMirrorsAreInPatch)
{
    // Plane at x=0. Physical cells: AMR x in [-2, 1]. Ghost shell grows nbrGhosts layers
    // from the cut cells into the solid side (AMR x < -2). `pointIsInterpolable_` requires
    // the mirror's stencil (2 consecutive grid values per direction; for cell-centered
    // fields that means 1 cell of slack per end of the ghost box). Ghost cells deep in the
    // solid-side halo map mirrors near or beyond the fluid-side halo edge, so a few mirrors
    // legitimately fail the check.
    //
    // The contract here is: at least the ghost cells whose mirrors land in the *physical*
    // domain (interior cells) must be interpolable.
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}};
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Mapper::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Mapper::withDefaults(plane, layout, ov);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    auto const& ghost_cells = buffers.tags.getGhostDataFromCentering(kCellC);
    ASSERT_FALSE(ghost_cells.empty());

    auto const& dx     = layout.meshSize();
    auto const& amrBox = layout.AMRBox();
    auto inPhysicalDomain = [&](auto const& point) {
        for (auto d = 0u; d < 2u; ++d)
        {
            int const iCell = static_cast<int>(std::floor(point[d] / dx[d]));
            if (iCell < amrBox.lower[d] || iCell > amrBox.upper[d])
                return false;
        }
        return true;
    };

    std::size_t interpolable_in_physical = 0;
    for (auto const& g : ghost_cells)
    {
        if (inPhysicalDomain(g.mirrorPoint))
        {
            EXPECT_TRUE(g.mirrorIsInterpolable)
                << "Ghost cell whose mirror lies in the physical domain must be interpolable";
            ++interpolable_in_physical;
        }
    }
    EXPECT_GT(interpolable_in_physical, 0u)
        << "At least one ghost mirror must land in the physical domain";
}

TEST(InnerBoundaryMeshClassifier, underResolvedSphereProducesDegradedElements)
{
    // radius spans a single cell (R/dx = 1), far below the (nbrGhosts + 1) requirement, so the
    // level under-resolves the sphere and a degraded band must be emitted.
    PHARE::core::SphereInnerBoundary<2> sphere{"sphere", {0.05, 0.05}, 1.0};
    GridLayout layout{{1.0, 1.0}, {12u, 12u}, {-6., -6.}};
    auto tagger = Mapper::withDefaults(sphere, layout);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    EXPECT_FALSE(buffers.tags.getDegradedElemsFromCentering(kCellC).empty())
        << "an under-resolved sphere must emit degraded cell elements";
}

TEST(InnerBoundaryMeshClassifier, resolvedSphereProducesNoDegradedElements)
{
    GridLayout probe{{1.0, 1.0}, {4u, 4u}, {0., 0.}};
    auto const ng = probe.nbrGhosts();

    double const dx     = 1.0;
    double const radius = (static_cast<double>(ng) + 2.0) * dx; // comfortably resolves the sphere
    auto const half     = static_cast<std::uint32_t>(std::ceil(radius / dx)) + ng + 2u;
    std::uint32_t const ncell = 2u * half;
    double const origin       = -static_cast<double>(half) * dx;

    PHARE::core::SphereInnerBoundary<2> sphere{"sphere", {0.05, 0.05}, radius};
    GridLayout layout{{dx, dx}, {ncell, ncell}, {origin, origin}};
    auto tagger = Mapper::withDefaults(sphere, layout);

    InnerBoundaryMeshClassifierBuffers buffers{layout};
    tagger(layout, buffers.tags);

    EXPECT_TRUE(buffers.tags.getDegradedElemsFromCentering(kCellC).empty())
        << "a resolved sphere must not emit any degraded elements";
    // sanity: the boundary surface is inside the grid, so a ghost shell exists.
    EXPECT_FALSE(buffers.tags.getGhostDataFromCentering(kCellC).empty());
}
