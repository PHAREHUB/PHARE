#include "gtest/gtest.h"

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/inner_boundary_mesh_mapper.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/inner_boundary/sphere_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"

namespace
{
constexpr double eps = 1e-12;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Mapper = PHARE::core::InnerBoundaryMeshMapper<2, GridLayout, PHARE::core::MHDQuantity>;
using MeshData           = PHARE::core::InnerBoundaryMeshData<2, PHARE::core::MHDQuantity>;
using NodeField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using CellField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using FaceField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using FaceVec        = PHARE::core::VecField<FaceField, PHARE::core::MHDQuantity>;
using EdgeField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using EdgeVec        = PHARE::core::VecField<EdgeField, PHARE::core::MHDQuantity>;

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

struct InnerBoundaryMeshMapperBuffers
{
    // InnerBoundaryMeshData names its fields as "<boundary>_<component>"; keep this in sync.
    static constexpr char const* BOUNDARY_NAME = "test";

    explicit InnerBoundaryMeshMapperBuffers(GridLayout const& layout)
        : phi_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::NodeCentered)}
        , cell_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)}
        , face_x_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::FaceCenteredX)}
        , face_y_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::FaceCenteredY)}
        , face_z_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::FaceCenteredZ)}
        , face_fields{FaceField{std::string(BOUNDARY_NAME) + "_face_status_x",
                                PHARE::core::MHDQuantity::Scalar::FaceCenteredX,
                                face_x_storage.data(), face_x_storage.shape()},
                      FaceField{std::string(BOUNDARY_NAME) + "_face_status_y",
                                PHARE::core::MHDQuantity::Scalar::FaceCenteredY,
                                face_y_storage.data(), face_y_storage.shape()},
                      FaceField{std::string(BOUNDARY_NAME) + "_face_status_z",
                                PHARE::core::MHDQuantity::Scalar::FaceCenteredZ,
                                face_z_storage.data(), face_z_storage.shape()}}
        , edge_x_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::EdgeCenteredX)}
        , edge_y_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::EdgeCenteredY)}
        , edge_z_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::EdgeCenteredZ)}
        , edge_fields{EdgeField{std::string(BOUNDARY_NAME) + "_edge_status_x",
                                PHARE::core::MHDQuantity::Scalar::EdgeCenteredX,
                                edge_x_storage.data(), edge_x_storage.shape()},
                      EdgeField{std::string(BOUNDARY_NAME) + "_edge_status_y",
                                PHARE::core::MHDQuantity::Scalar::EdgeCenteredY,
                                edge_y_storage.data(), edge_y_storage.shape()},
                      EdgeField{std::string(BOUNDARY_NAME) + "_edge_status_z",
                                PHARE::core::MHDQuantity::Scalar::EdgeCenteredZ,
                                edge_z_storage.data(), edge_z_storage.shape()}}
        , tags{BOUNDARY_NAME}
    {
        NodeField phi_field{std::string(BOUNDARY_NAME) + "_signed_distance",
                            PHARE::core::MHDQuantity::Scalar::NodeCentered,
                            phi_storage.data(), phi_storage.shape()};
        tags.signedDistanceAtNodes.setBuffer(&phi_field);

        CellField cell_field{std::string(BOUNDARY_NAME) + "_cell_status",
                             PHARE::core::MHDQuantity::Scalar::CellCentered,
                             cell_storage.data(), cell_storage.shape()};
        tags.cellStatus.setBuffer(&cell_field);

        tags.faceStatus.setBuffer(&face_fields);
        tags.edgeStatus.setBuffer(&edge_fields);
    }

    PHARE::core::NdArrayVector<2, double> phi_storage;

    PHARE::core::NdArrayVector<2, double> cell_storage;

    PHARE::core::NdArrayVector<2, double> face_x_storage;
    PHARE::core::NdArrayVector<2, double> face_y_storage;
    PHARE::core::NdArrayVector<2, double> face_z_storage;
    std::array<FaceField, 3> face_fields;

    PHARE::core::NdArrayVector<2, double> edge_x_storage;
    PHARE::core::NdArrayVector<2, double> edge_y_storage;
    PHARE::core::NdArrayVector<2, double> edge_z_storage;
    std::array<EdgeField, 3> edge_fields;
    MeshData tags;
};
} // namespace

TEST(InnerBoundaryMeshMapper, computesReasonableDefaultCutEpsFromLayout)
{
    PHARE::core::SphereInnerBoundary<2> sphere{"sphere", {0., 0.}, 1.};
    GridLayout layout{{0.2, 0.1}, {4u, 4u}, {0., 0.}};
    auto tagger = Mapper::withDefaults(sphere, layout);

    InnerBoundaryMeshMapperBuffers buffers{layout};
    tagger(layout, buffers.tags);

    auto const origin_node = physicalLocalIndex(layout, buffers.tags.signedDistanceAtNodes, 0u, 0u);
    EXPECT_NEAR(buffers.tags.signedDistanceAtNodes(origin_node), -1., eps);
}

TEST(InnerBoundaryMeshMapper, tagsCutInactiveAndGhostGeometry)
{
    PHARE::core::PlaneInnerBoundary<2> plane{"plane", {0.0, 0.0}, {1.0, 0.0}}; // x=0
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Mapper::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Mapper::withDefaults(plane, layout, ov);

    InnerBoundaryMeshMapperBuffers buffers{layout};
    tagger(layout, buffers.tags);

    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::CellStatus::Inactive),
              buffers.tags.cellStatus(physicalLocalIndex(layout, buffers.tags.cellStatus, 0u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::CellStatus::Cut),
              buffers.tags.cellStatus(physicalLocalIndex(layout, buffers.tags.cellStatus, 1u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::CellStatus::Cut),
              buffers.tags.cellStatus(physicalLocalIndex(layout, buffers.tags.cellStatus, 2u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::CellStatus::Ghost),
              buffers.tags.cellStatus(physicalLocalIndex(layout, buffers.tags.cellStatus, 3u, 0u)));

    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::FaceStatus::Cut),
              buffers.tags.faceStatus[0](physicalLocalIndex(layout, buffers.tags.faceStatus[0], 2u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::FaceStatus::Ghost),
              buffers.tags.faceStatus[0](physicalLocalIndex(layout, buffers.tags.faceStatus[0], 3u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::FaceStatus::Ghost),
              buffers.tags.faceStatus[1](physicalLocalIndex(layout, buffers.tags.faceStatus[1], 3u, 0u)));

    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::EdgeStatus::Cut),
              buffers.tags.edgeStatus[1](physicalLocalIndex(layout, buffers.tags.edgeStatus[1], 2u, 0u)));
    EXPECT_EQ(PHARE::core::toDouble(PHARE::core::EdgeStatus::Ghost),
              buffers.tags.edgeStatus[1](physicalLocalIndex(layout, buffers.tags.edgeStatus[1], 3u, 0u)));
}
