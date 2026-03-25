#include "gtest/gtest.h"

#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee_mhd.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/inner_boundary/inner_boundary_geometry_tagger.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/inner_boundary/sphere_inner_boundary.hpp"
#include "core/mhd/mhd_quantities.hpp"

namespace
{
constexpr double eps = 1e-12;

using GridLayoutImpl = PHARE::core::GridLayoutImplYeeMHD<2, 2, 1>;
using GridLayout     = PHARE::core::GridLayout<GridLayoutImpl>;
using Tagger = PHARE::core::InnerBoundaryGeometryTagger<2, GridLayout, PHARE::core::MHDQuantity>;
using Tags           = PHARE::core::InnerBoundaryTags<2, PHARE::core::MHDQuantity>;
using NodeField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, double>;
using CellField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, PHARE::core::CellTag>;
using FaceField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, PHARE::core::FaceTag>;
using FaceVec        = PHARE::core::VecField<FaceField, PHARE::core::MHDQuantity>;
using EdgeField      = PHARE::core::Field<2, PHARE::core::MHDQuantity::Scalar, PHARE::core::EdgeTag>;
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

struct InnerBoundaryGeometryTaggerBuffers
{
    explicit InnerBoundaryGeometryTaggerBuffers(GridLayout const& layout)
        : phi_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::NodeCentered)}
        , cell_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::CellCentered)}
        , face_x_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::FaceCenteredX)}
        , face_y_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::FaceCenteredY)}
        , face_z_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::FaceCenteredZ)}
        , face_fields{
              FaceField{"face_tags_x", PHARE::core::MHDQuantity::Scalar::FaceCenteredX,
                        face_x_storage.data(), face_x_storage.shape()},
               FaceField{"face_tags_y", PHARE::core::MHDQuantity::Scalar::FaceCenteredY,
                         face_y_storage.data(), face_y_storage.shape()},
               FaceField{"face_tags_z", PHARE::core::MHDQuantity::Scalar::FaceCenteredZ,
                         face_z_storage.data(), face_z_storage.shape()}}
        , edge_x_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::EdgeCenteredX)}
        , edge_y_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::EdgeCenteredY)}
        , edge_z_storage{layout.allocSize(PHARE::core::MHDQuantity::Scalar::EdgeCenteredZ)}
        , edge_fields{
              EdgeField{"edge_tags_x", PHARE::core::MHDQuantity::Scalar::EdgeCenteredX,
                        edge_x_storage.data(), edge_x_storage.shape()},
              EdgeField{"edge_tags_y", PHARE::core::MHDQuantity::Scalar::EdgeCenteredY,
                        edge_y_storage.data(), edge_y_storage.shape()},
              EdgeField{"edge_tags_z", PHARE::core::MHDQuantity::Scalar::EdgeCenteredZ,
                        edge_z_storage.data(), edge_z_storage.shape()}}
        , tags{NodeField{"phi_nodes", PHARE::core::MHDQuantity::Scalar::NodeCentered,
                         phi_storage.data(), phi_storage.shape()},
               CellField{"cell_tags", PHARE::core::MHDQuantity::Scalar::CellCentered,
                         cell_storage.data(), cell_storage.shape()},
               FaceVec{"face_tags", PHARE::core::MHDQuantity::Vector::FaceCentered},
               EdgeVec{"edge_tags", PHARE::core::MHDQuantity::Vector::EdgeCentered}}
    {
        tags.faceTags.setBuffer(&face_fields);
        tags.edgeTags.setBuffer(&edge_fields);
    }

    PHARE::core::NdArrayVector<2, double> phi_storage;

    PHARE::core::NdArrayVector<2, PHARE::core::CellTag> cell_storage;

    PHARE::core::NdArrayVector<2, PHARE::core::FaceTag> face_x_storage;
    PHARE::core::NdArrayVector<2, PHARE::core::FaceTag> face_y_storage;
    PHARE::core::NdArrayVector<2, PHARE::core::FaceTag> face_z_storage;
    std::array<FaceField, 3> face_fields;

    PHARE::core::NdArrayVector<2, PHARE::core::EdgeTag> edge_x_storage;
    PHARE::core::NdArrayVector<2, PHARE::core::EdgeTag> edge_y_storage;
    PHARE::core::NdArrayVector<2, PHARE::core::EdgeTag> edge_z_storage;
    std::array<EdgeField, 3> edge_fields;
    Tags tags;
};
} // namespace

TEST(InnerBoundaryGeometryTagger, computesReasonableDefaultCutEpsFromLayout)
{
    PHARE::core::SphereInnerBoundary<2> sphere{{0., 0.}, 1.};
    GridLayout layout{{0.2, 0.1}, {4u, 4u}, {0., 0.}};
    auto tagger = Tagger::withDefaults(sphere, layout);

    InnerBoundaryGeometryTaggerBuffers buffers{layout};
    tagger.tagAll(layout, buffers.tags);

    auto const origin_node = physicalLocalIndex(layout, buffers.tags.phiNodes, 0u, 0u);
    EXPECT_NEAR(buffers.tags.phiNodes(origin_node), -1., eps);
}

TEST(InnerBoundaryGeometryTagger, tagsCutInactiveAndGhostGeometry)
{
    PHARE::core::PlaneInnerBoundary<2> plane{{0.0, 0.0}, {1.0, 0.0}}; // x=0
    PHARE::core::Box<int, 2> amr_box{{-2, 0}, {1, 1}};
    GridLayout layout{{1.0, 1.0}, {4u, 2u}, {0.0, 0.0}, amr_box};

    Tagger::Overrides ov;
    ov.cut_eps      = 1e-12;
    ov.inactive_eps = 0.0;
    auto tagger     = Tagger::withDefaults(plane, layout, ov);

    InnerBoundaryGeometryTaggerBuffers buffers{layout};
    tagger.tagAll(layout, buffers.tags);

    EXPECT_EQ(PHARE::core::CellTag::Inactive,
              buffers.tags.cellTags(physicalLocalIndex(layout, buffers.tags.cellTags, 0u, 0u)));
    EXPECT_EQ(PHARE::core::CellTag::Cut,
              buffers.tags.cellTags(physicalLocalIndex(layout, buffers.tags.cellTags, 1u, 0u)));
    EXPECT_EQ(PHARE::core::CellTag::Cut,
              buffers.tags.cellTags(physicalLocalIndex(layout, buffers.tags.cellTags, 2u, 0u)));
    EXPECT_EQ(PHARE::core::CellTag::Ghost,
              buffers.tags.cellTags(physicalLocalIndex(layout, buffers.tags.cellTags, 3u, 0u)));

    EXPECT_EQ(PHARE::core::FaceTag::Cut,
              buffers.tags.faceTags[0](physicalLocalIndex(layout, buffers.tags.faceTags[0], 2u, 0u)));
    EXPECT_EQ(PHARE::core::FaceTag::Ghost,
              buffers.tags.faceTags[0](physicalLocalIndex(layout, buffers.tags.faceTags[0], 3u, 0u)));
    EXPECT_EQ(PHARE::core::FaceTag::Ghost,
              buffers.tags.faceTags[1](physicalLocalIndex(layout, buffers.tags.faceTags[1], 3u, 0u)));

    EXPECT_EQ(PHARE::core::EdgeTag::Cut,
              buffers.tags.edgeTags[1](physicalLocalIndex(layout, buffers.tags.edgeTags[1], 2u, 0u)));
    EXPECT_EQ(PHARE::core::EdgeTag::Ghost,
              buffers.tags.edgeTags[1](physicalLocalIndex(layout, buffers.tags.edgeTags[1], 3u, 0u)));
}
