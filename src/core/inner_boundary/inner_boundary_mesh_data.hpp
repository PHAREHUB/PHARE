#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_DATA_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_DATA_HPP

#include "core/data/field/field.hpp"
#include "core/data/vecfield/vecfield.hpp"

#include <cstdint>

namespace PHARE::core
{
/**
 * @brief Cell status relative to an inner boundary.
 */
enum class CellStatus : std::uint8_t { Fluid, Cut, Ghost, Inactive };
/**
 * @brief Face status relative to an inner boundary.
 */
enum class FaceStatus : std::uint8_t { Fluid, Cut, Ghost, Inactive };
/**
 * @brief Edge status relative to an inner boundary.
 */
enum class EdgeStatus : std::uint8_t { Fluid, Cut, Ghost, Inactive };

/// Convert a CellStatus enum value to its double encoding for field storage.
inline constexpr double toDouble(CellStatus s)
{
    return static_cast<double>(static_cast<std::uint8_t>(s));
}
/// Convert a FaceStatus enum value to its double encoding for field storage.
inline constexpr double toDouble(FaceStatus s)
{
    return static_cast<double>(static_cast<std::uint8_t>(s));
}
/// Convert an EdgeStatus enum value to its double encoding for field storage.
inline constexpr double toDouble(EdgeStatus s)
{
    return static_cast<double>(static_cast<std::uint8_t>(s));
}


/**
 * @brief Gather useful type definitions for mesh data around an inner boundary.
 */
template<std::size_t dim, typename PhysicalQuantityT>
struct InnerBoundaryMeshDataTypes
{
    using node_double_field_type    = Field<dim, typename PhysicalQuantityT::Scalar, double>;
    using cell_status_field_type    = Field<dim, typename PhysicalQuantityT::Scalar, double>;
    using face_status_field_type    = Field<dim, typename PhysicalQuantityT::Scalar, double>;
    using face_status_vecfield_type = VecField<face_status_field_type, PhysicalQuantityT>;
    using edge_status_field_type    = Field<dim, typename PhysicalQuantityT::Scalar, double>;
    using edge_status_vecfield_type = VecField<edge_status_field_type, PhysicalQuantityT>;
};

/**
 * @brief Bundle of node level-set values and mesh status data around an inner boundary.
 *
 * @tparam dim Spatial dimension.
 * @tparam PhysicalQuantityT Quantity traits used to define node/cell/face/edge
 * field types.
 */
template<std::size_t dim, typename PhysicalQuantityT>
struct InnerBoundaryMeshData
{
    using types                     = InnerBoundaryMeshDataTypes<dim, PhysicalQuantityT>;
    using node_double_field_type    = types::node_double_field_type;
    using cell_status_field_type    = types::cell_status_field_type;
    using face_status_vecfield_type = types::face_status_vecfield_type;
    using edge_status_vecfield_type = types::edge_status_vecfield_type;


    InnerBoundaryMeshData(std::string boundaryName = "")
        : signedDistanceAtNodes{boundaryName + "_signed_distance",
                                PhysicalQuantityT::Scalar::NodeCentered}
        , cellStatus{boundaryName + "_cell_status", PhysicalQuantityT::Scalar::CellCentered}
        , faceStatus{boundaryName + "_face_status", PhysicalQuantityT::Vector::FaceCentered}
        , edgeStatus{boundaryName + "_edge_status", PhysicalQuantityT::Vector::EdgeCentered}
    {
    }

    node_double_field_type signedDistanceAtNodes; ///< Signed distance to the boundary at nodes.
    cell_status_field_type cellStatus;            ///< Per-cell fluid/cut/ghost/inactive
                                                  ///< classification.
    face_status_vecfield_type faceStatus;         ///< Per-face fluid/cut/ghost/inactive
                                                  ///< classification.
    edge_status_vecfield_type edgeStatus;         ///< Per-edge fluid/cut/ghost/inactive
                                                  ///< classification.



    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return signedDistanceAtNodes.isUsable() and cellStatus.isUsable() and faceStatus.isUsable()
               and edgeStatus.isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return signedDistanceAtNodes.isSettable() and cellStatus.isSettable()
               and faceStatus.isSettable() and edgeStatus.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(signedDistanceAtNodes, cellStatus, faceStatus, edgeStatus);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(signedDistanceAtNodes, cellStatus, faceStatus, edgeStatus);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------
};
} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_DATA_HPP
