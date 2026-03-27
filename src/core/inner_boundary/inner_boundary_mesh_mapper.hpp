#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_MAPPER_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_MAPPER_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/inner_boundary/inner_boundary.hpp"
#include "core/inner_boundary/inner_boundary_mesh_data.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace PHARE::core
{



/**
 * @brief Tags cells, faces, and edges around an embedded inner boundary.
 *
 * The tagger first evaluates the boundary signed distance on the node support,
 * classifies cells as fluid/cut/inactive, then grows a ghost-cell shell around
 * cut cells and finally propagates that geometry to faces and edges.
 *
 * @tparam dim Spatial dimension.
 * @tparam GridLayoutT Grid layout type used to provide coordinates and iterate
 * local supports.
 * @tparam PhysicalQuantityT Quantity traits used to define node/cell/face/edge
 * field types.
 */
template<std::size_t dim, typename GridLayoutT, typename PhysicalQuantityT>
class InnerBoundaryMeshMapper
{
public:
    using mesh_data_types           = InnerBoundaryMeshDataTypes<dim, PhysicalQuantityT>;
    using point_type                = typename InnerBoundary<dim>::point_type;
    using local_index_type          = Point<std::uint32_t, dim>;
    using signed_local_index_type   = Point<int, dim>;
    using node_double_field_type    = mesh_data_types::node_double_field_type;
    using cell_status_field_type    = mesh_data_types::cell_status_field_type;
    using face_status_field_type    = mesh_data_types::face_status_field_type;
    using face_status_vecfield_type = mesh_data_types::face_status_vecfield_type;
    using edge_status_field_type    = mesh_data_types::edge_status_field_type;
    using edge_status_vecfield_type = mesh_data_types::edge_status_vecfield_type;
    using mesh_data_type            = InnerBoundaryMeshData<dim, PhysicalQuantityT>;


    /**
     * @brief Runtime parameters controlling cut and ghost classification.
     *
     * @param ghost_layers Number of fluid-cell layers promoted to ghost around
     * cut cells.
     * @param cut_eps Tolerance used when deciding whether a support intersects
     * the boundary.
     * @param inactive_eps Tolerance used when deciding whether a non-cut support
     * lies inside the boundary.
     */
    struct Params
    {
        std::size_t ghost_layers = 1;
        double cut_eps           = 0.0;
        double inactive_eps      = 0.0;
    };

    /**
     * @brief Optional overrides for default parameter construction.
     *
     * Any unset entry keeps the mesh-derived default.
     */
    struct Overrides
    {
        std::optional<std::size_t> ghost_layers;
        std::optional<double> cut_eps;
        std::optional<double> inactive_eps;
    };

    /**
     * @brief Build a tagger with mesh-derived default tolerances.
     *
     * @param boundary Embedded boundary used for signed-distance queries.
     * @param layout Grid layout used to derive default tolerances from the mesh
     * spacing.
     * @param overrides Optional parameter overrides.
     * @return Configured tagger instance.
     */
    static InnerBoundaryMeshMapper withDefaults(InnerBoundary<dim> const& boundary,
                                                GridLayoutT const& layout,
                                                Overrides const& overrides = {})
    {
        auto const& dx    = layout.meshSize();
        auto const dx_min = *std::min_element(dx.begin(), dx.end());

        Params p;
        p.ghost_layers = overrides.ghost_layers.value_or(1);
        p.cut_eps      = overrides.cut_eps.value_or(1e-6 * dx_min);
        p.inactive_eps = overrides.inactive_eps.value_or(p.cut_eps);
        return InnerBoundaryMeshMapper{boundary, p};
    }

    /**
     * @brief Construct a tagger with explicit parameters.
     *
     * @param boundary Embedded boundary used for signed-distance queries.
     * @param params Explicit cut/ghost classification parameters.
     */
    InnerBoundaryMeshMapper(InnerBoundary<dim> const& boundary, Params params)
        : boundary_{boundary}
        , params_{params}
    {
    }

    /**
     * @brief Fill node level-set values and classify cells, faces, and edges.
     *
     * @param layout Grid layout used to iterate and locate all supports.
     * @param signed_distance_at_nodes Node-centered signed-distance field to populate.
     * @param cell_status Cell meshData written by the tagger.
     * @param face_status Face meshData written by the tagger.
     * @param edge_status Edge meshData written by the tagger.
     */
    void operator()(GridLayoutT const& layout, node_double_field_type& signed_distance_at_nodes,
                    cell_status_field_type& cell_status, face_status_vecfield_type& face_status,
                    edge_status_vecfield_type& edge_status) const
    {
        validateCenterings_(signed_distance_at_nodes, cell_status, face_status, edge_status);
        fillNodePhi_(layout, signed_distance_at_nodes);
        tagCutInactiveAndFluidCells_(layout, signed_distance_at_nodes, cell_status);
        tagGhostCells_(layout, cell_status);
        tagFaces_(layout, signed_distance_at_nodes, cell_status, face_status);
        tagEdges_(layout, signed_distance_at_nodes, cell_status, edge_status);
    }

    /**
     * @brief Fill node level-set values and classify all supports in a tag bundle.
     *
     * @param layout Grid layout used to iterate and locate all supports.
     * @param meshData Bundle written by the tagger.
     */
    void operator()(GridLayoutT const& layout, mesh_data_type& meshData) const
    {
        (*this)(layout, meshData.signedDistanceAtNodes, meshData.cellStatus, meshData.faceStatus,
                meshData.edgeStatus);
    }

private:
    InnerBoundary<dim> const& boundary_;
    Params params_;

    /**
     * @brief Decide whether a support is geometrically cut by the boundary.
     *
     * @param phi_min Minimum signed distance on the sampled support.
     * @param phi_max Maximum signed distance on the sampled support.
     * @param cut_eps Cut tolerance.
     * @return True when the sampled support intersects the boundary.
     */
    static bool isCut_(double phi_min, double phi_max, double cut_eps)
    {
        if (phi_min < -cut_eps && phi_max > cut_eps)
            return true;
        return phi_min <= cut_eps && phi_max >= -cut_eps;
    }

    /**
     * @brief Check whether a signed local index lies inside a local array shape.
     *
     * @param idx Candidate signed local index.
     * @param shape Array extents.
     * @return True when @p idx is a valid local index for @p shape.
     */
    static bool inBounds_(signed_local_index_type const& idx,
                          std::array<std::uint32_t, dim> const& shape)
    {
        for (std::size_t i = 0; i < dim; ++i)
            if (idx[i] < 0 || idx[i] >= static_cast<int>(shape[i]))
                return false;
        return true;
    }

    /**
     * @brief Convert an unsigned local index to its signed counterpart.
     *
     * @param idx Unsigned local index.
     * @return Signed version of @p idx.
     */
    static signed_local_index_type asSigned_(local_index_type const& idx)
    {
        signed_local_index_type signed_idx{};
        for (std::size_t i = 0; i < dim; ++i)
            signed_idx[i] = static_cast<int>(idx[i]);
        return signed_idx;
    }

    /**
     * @brief Convert a signed local index back to an unsigned local index.
     *
     * Caller code is expected to have checked bounds beforehand.
     *
     * @param idx Signed local index.
     * @return Unsigned version of @p idx.
     */
    static local_index_type asLocal_(signed_local_index_type const& idx)
    {
        local_index_type local_idx{};
        for (std::size_t i = 0; i < dim; ++i)
            local_idx[i] = static_cast<std::uint32_t>(idx[i]);
        return local_idx;
    }

    /**
     * @brief Map a field-local support index to the corresponding AMR index.
     *
     * @tparam FieldT Field type whose centering determines the local-to-AMR
     * offset.
     * @param layout Grid layout providing AMR and physical-start information.
     * @param field Field whose centering is being mapped.
     * @param local_idx Local field index.
     * @return AMR index matching @p local_idx for @p field.
     */
    template<typename FieldT>
    static signed_local_index_type fieldAMRIndex_(GridLayoutT const& layout, FieldT const& field,
                                                  local_index_type const& local_idx)
    {
        // GridLayout::localToAMR is defined for cell-centered indices only. The
        // tagger also manipulates node-, face-, and edge-centered supports, so
        // the shift must be recomputed from the field centering.
        signed_local_index_type amr_idx{};
        for (std::size_t i = 0; i < dim; ++i)
        {
            auto const dir = static_cast<Direction>(i);
            amr_idx[i]     = static_cast<int>(local_idx[i])
                         + (layout.AMRBox().lower[i]
                            - static_cast<int>(layout.physicalStartIndex(field, dir)));
        }
        return amr_idx;
    }

    /**
     * @brief Expected centering for node-based supports.
     *
     * @return Node/primal centering descriptor.
     */
    static auto nodeCentering_()
    {
        std::array<QtyCentering, dim> centering{};
        centering.fill(QtyCentering::primal);
        return centering;
    }

    /**
     * @brief Expected centering for cell-based supports.
     *
     * @return Cell/dual centering descriptor.
     */
    static auto cellCentering_()
    {
        std::array<QtyCentering, dim> centering{};
        centering.fill(QtyCentering::dual);
        return centering;
    }

    /**
     * @brief Expected centering for the face family aligned with one direction.
     *
     * @param dir Face normal direction.
     * @return Centering descriptor for faces normal to @p dir.
     */
    static auto faceCentering_(std::size_t dir)
    {
        auto centering = cellCentering_();
        centering[dir] = QtyCentering::primal;
        return centering;
    }

    /**
     * @brief Expected centering for the edge family aligned with one direction.
     *
     * @param dir Edge direction.
     * @return Centering descriptor for edges aligned with @p dir.
     */
    static auto edgeCentering_(std::size_t dir)
    {
        auto centering = nodeCentering_();
        centering[dir] = QtyCentering::dual;
        return centering;
    }

    /**
     * @brief Visit the node corners of a cell-local support.
     *
     * @tparam iDim Current recursion dimension.
     * @tparam Fn Callable receiving each corner index.
     * @param cell Cell local index.
     * @param node Scratch node index mutated during recursion.
     * @param fn Callback invoked for each corner.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachCellCorner_(local_index_type const& cell, local_index_type& node, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(node);
        }
        else
        {
            node[iDim] = cell[iDim];
            forEachCellCorner_<iDim + 1>(cell, node, std::forward<Fn>(fn));
            node[iDim] = cell[iDim] + 1;
            forEachCellCorner_<iDim + 1>(cell, node, std::forward<Fn>(fn));
        }
    }

    /**
     * @brief Visit the node corners of a face-local support.
     *
     * @tparam iDim Current recursion dimension.
     * @tparam Fn Callable receiving each corner index.
     * @param face Face local index.
     * @param dir Face normal direction.
     * @param node Scratch node index mutated during recursion.
     * @param fn Callback invoked for each corner.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachFaceCorner_(local_index_type const& face, std::size_t dir,
                                   local_index_type& node, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(node);
        }
        else if (iDim == dir)
        {
            node[iDim] = face[iDim];
            forEachFaceCorner_<iDim + 1>(face, dir, node, std::forward<Fn>(fn));
        }
        else
        {
            node[iDim] = face[iDim];
            forEachFaceCorner_<iDim + 1>(face, dir, node, std::forward<Fn>(fn));
            node[iDim] = face[iDim] + 1;
            forEachFaceCorner_<iDim + 1>(face, dir, node, std::forward<Fn>(fn));
        }
    }

    /**
     * @brief Visit the node endpoints of an edge-local support.
     *
     * @tparam iDim Current recursion dimension.
     * @tparam Fn Callable receiving each endpoint index.
     * @param edge Edge local index.
     * @param dir Edge direction.
     * @param node Scratch node index mutated during recursion.
     * @param fn Callback invoked for each endpoint.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachEdgeEndpoint_(local_index_type const& edge, std::size_t dir,
                                     local_index_type& node, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(node);
        }
        else if (iDim == dir)
        {
            node[iDim] = edge[iDim];
            forEachEdgeEndpoint_<iDim + 1>(edge, dir, node, std::forward<Fn>(fn));
            node[iDim] = edge[iDim] + 1;
            forEachEdgeEndpoint_<iDim + 1>(edge, dir, node, std::forward<Fn>(fn));
        }
        else
        {
            node[iDim] = edge[iDim];
            forEachEdgeEndpoint_<iDim + 1>(edge, dir, node, std::forward<Fn>(fn));
        }
    }

    /**
     * @brief Visit the cells adjacent to a face-local support.
     *
     * @tparam iDim Current recursion dimension.
     * @tparam Fn Callable receiving each adjacent cell index.
     * @param face Face local index.
     * @param dir Face normal direction.
     * @param cell Scratch signed cell index mutated during recursion.
     * @param fn Callback invoked for each adjacent cell.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachFaceAdjacentCell_(local_index_type const& face, std::size_t dir,
                                         signed_local_index_type& cell, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(cell);
        }
        else if (iDim == dir)
        {
            cell[iDim] = static_cast<int>(face[iDim]) - 1;
            forEachFaceAdjacentCell_<iDim + 1>(face, dir, cell, std::forward<Fn>(fn));
            cell[iDim] = static_cast<int>(face[iDim]);
            forEachFaceAdjacentCell_<iDim + 1>(face, dir, cell, std::forward<Fn>(fn));
        }
        else
        {
            cell[iDim] = static_cast<int>(face[iDim]);
            forEachFaceAdjacentCell_<iDim + 1>(face, dir, cell, std::forward<Fn>(fn));
        }
    }

    /**
     * @brief Visit the cells adjacent to an edge-local support.
     *
     * @tparam iDim Current recursion dimension.
     * @tparam Fn Callable receiving each adjacent cell index.
     * @param edge Edge local index.
     * @param dir Edge direction.
     * @param cell Scratch signed cell index mutated during recursion.
     * @param fn Callback invoked for each adjacent cell.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachEdgeAdjacentCell_(local_index_type const& edge, std::size_t dir,
                                         signed_local_index_type& cell, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(cell);
        }
        else if (iDim == dir)
        {
            cell[iDim] = static_cast<int>(edge[iDim]);
            forEachEdgeAdjacentCell_<iDim + 1>(edge, dir, cell, std::forward<Fn>(fn));
        }
        else
        {
            cell[iDim] = static_cast<int>(edge[iDim]) - 1;
            forEachEdgeAdjacentCell_<iDim + 1>(edge, dir, cell, std::forward<Fn>(fn));
            cell[iDim] = static_cast<int>(edge[iDim]);
            forEachEdgeAdjacentCell_<iDim + 1>(edge, dir, cell, std::forward<Fn>(fn));
        }
    }

    /**
     * @brief Check whether a face borders at least one ghost cell.
     *
     * @param face Face local index.
     * @param dir Face normal direction.
     * @param cell_status Cell meshData used to inspect adjacent cells.
     * @return True when the face surrounds a ghost cell.
     */
    bool faceSurroundsGhostCell_(local_index_type const& face, std::size_t dir,
                                 cell_status_field_type const& cell_status) const
    {
        bool has_ghost        = false;
        auto const cell_shape = cell_status.shape();
        signed_local_index_type cell{};
        forEachFaceAdjacentCell_<0>(face, dir, cell, [&](auto const& adjacent_cell) {
            if (has_ghost || !inBounds_(adjacent_cell, cell_shape))
                return;

            if (cell_status(asLocal_(adjacent_cell)) == toDouble(CellStatus::Ghost))
                has_ghost = true;
        });
        return has_ghost;
    }

    /**
     * @brief Check whether an edge borders at least one ghost cell.
     *
     * @param edge Edge local index.
     * @param dir Edge direction.
     * @param cell_status Cell meshData used to inspect adjacent cells.
     * @return True when the edge surrounds a ghost cell.
     */
    bool edgeSurroundsGhostCell_(local_index_type const& edge, std::size_t dir,
                                 cell_status_field_type const& cell_status) const
    {
        bool has_ghost        = false;
        auto const cell_shape = cell_status.shape();
        signed_local_index_type cell{};
        forEachEdgeAdjacentCell_<0>(edge, dir, cell, [&](auto const& adjacent_cell) {
            if (has_ghost || !inBounds_(adjacent_cell, cell_shape))
                return;

            if (cell_status(asLocal_(adjacent_cell)) == toDouble(CellStatus::Ghost))
                has_ghost = true;
        });
        return has_ghost;
    }

    /**
     * @brief Populate the node-centered signed-distance field.
     *
     * @param layout Grid layout used to iterate and locate nodes.
     * @param signed_distance_at_nodes Node-centered level-set field to fill.
     */
    void fillNodePhi_(GridLayoutT const& layout,
                      node_double_field_type& signed_distance_at_nodes) const
    {
        layout.evalOnGhostBox(signed_distance_at_nodes, [&](auto... idx) {
            auto const local_point = local_index_type{static_cast<std::uint32_t>(idx)...};
            auto const amr_point   = fieldAMRIndex_(layout, signed_distance_at_nodes, local_point);
            signed_distance_at_nodes(local_point) = boundary_.signedDistance(
                layout.fieldNodeCoordinates(signed_distance_at_nodes, amr_point));
        });
    }

    /**
     * @brief Grow ghost-cell layers outward from cut cells.
     *
     * @param layout Grid layout used to iterate the local cell support.
     * @param cell_status Cell meshData updated in place.
     */
    void tagGhostCells_(GridLayoutT const& layout, cell_status_field_type& cell_status) const
    {
        if (params_.ghost_layers == 0)
            return;

        auto const cell_shape = cell_status.shape();
        std::vector<signed_local_index_type> frontier;

        layout.evalOnGhostBox(cell_status, [&](auto... idx) {
            auto const cell = local_index_type{static_cast<std::uint32_t>(idx)...};
            if (cell_status(cell) == toDouble(CellStatus::Cut))
                frontier.push_back(asSigned_(cell));
        });

        for (std::size_t layer = 0; layer < params_.ghost_layers && !frontier.empty(); ++layer)
        {
            std::vector<signed_local_index_type> next_frontier;
            next_frontier.reserve(frontier.size() * 2 * dim);

            // Grow a Manhattan shell outward from the previous frontier while
            // leaving cut and inactive cells untouched.
            for (auto const& source : frontier)
            {
                for (std::size_t d = 0; d < dim; ++d)
                {
                    for (int sgn : {-1, 1})
                    {
                        auto neigh = source;
                        neigh[d] += sgn;

                        if (!inBounds_(neigh, cell_shape))
                            continue;

                        auto const local_neigh = asLocal_(neigh);
                        if (cell_status(local_neigh) != toDouble(CellStatus::Fluid))
                            continue;

                        cell_status(local_neigh) = toDouble(CellStatus::Ghost);
                        next_frontier.push_back(neigh);
                    }
                }
            }

            frontier = std::move(next_frontier);
        }
    }

    /**
     * @brief Classify cells as cut, inactive, or fluid from geometry only.
     *
     * This pass intentionally does not assign ghost meshData.
     *
     * @param layout Grid layout used to iterate and locate cells.
     * @param signed_distance_at_nodes Node-centered signed-distance field.
     * @param cell_status Cell meshData updated in place.
     */
    void tagCutInactiveAndFluidCells_(GridLayoutT const& layout,
                                      node_double_field_type const& signed_distance_at_nodes,
                                      cell_status_field_type& cell_status) const
    {
        layout.evalOnGhostBox(cell_status, [&](auto... idx) {
            auto const local_cell = local_index_type{static_cast<std::uint32_t>(idx)...};
            auto const amr_cell   = fieldAMRIndex_(layout, cell_status, local_cell);

            double phi_min = std::numeric_limits<double>::max();
            double phi_max = std::numeric_limits<double>::lowest();
            local_index_type node_idx{};
            forEachCellCorner_<0>(local_cell, node_idx, [&](auto const& corner) {
                auto const phi = signed_distance_at_nodes(corner);
                phi_min        = std::min(phi_min, phi);
                phi_max        = std::max(phi_max, phi);
            });

            if (isCut_(phi_min, phi_max, params_.cut_eps))
            {
                cell_status(local_cell) = toDouble(CellStatus::Cut);
                return;
            }

            auto const phi_center
                = boundary_.signedDistance(layout.cellCenteredCoordinates(amr_cell));
            cell_status(local_cell)
                = (phi_center < -params_.inactive_eps) ? toDouble(CellStatus::Inactive) : toDouble(CellStatus::Fluid);
        });
    }

    /**
     * @brief Classify faces from geometry and cell ghost information.
     *
     * @param layout Grid layout used to iterate and locate faces.
     * @param signed_distance_at_nodes Node-centered signed-distance field.
     * @param cell_status Cell meshData used to detect ghost adjacency.
     * @param face_status Face meshData updated in place.
     */
    void tagFaces_(GridLayoutT const& layout,
                   node_double_field_type const& signed_distance_at_nodes,
                   cell_status_field_type const& cell_status,
                   face_status_vecfield_type& face_status) const
    {
        for (std::size_t dir = 0; dir < dim; ++dir)
        {
            layout.evalOnGhostBox(face_status[dir], [&](auto... idx) {
                auto const local_face = local_index_type{static_cast<std::uint32_t>(idx)...};
                auto const amr_face   = fieldAMRIndex_(layout, face_status[dir], local_face);

                double phi_min = std::numeric_limits<double>::max();
                double phi_max = std::numeric_limits<double>::lowest();
                local_index_type node_idx{};
                forEachFaceCorner_<0>(local_face, dir, node_idx, [&](auto const& corner) {
                    auto const phi = signed_distance_at_nodes(corner);
                    phi_min        = std::min(phi_min, phi);
                    phi_max        = std::max(phi_max, phi);
                });

                if (isCut_(phi_min, phi_max, params_.cut_eps))
                {
                    face_status[dir](local_face) = toDouble(FaceStatus::Cut);
                    return;
                }

                // Cut faces take priority; otherwise, any face adjacent to the
                // ghost-cell shell is itself tagged as ghost.
                if (faceSurroundsGhostCell_(local_face, dir, cell_status))
                {
                    face_status[dir](local_face) = toDouble(FaceStatus::Ghost);
                    return;
                }

                auto const phi_fc = boundary_.signedDistance(
                    layout.fieldNodeCoordinates(face_status[dir], amr_face));
                face_status[dir](local_face)
                    = (phi_fc < -params_.inactive_eps) ? toDouble(FaceStatus::Inactive) : toDouble(FaceStatus::Fluid);
            });
        }
    }

    /**
     * @brief Classify edges from geometry and cell ghost information.
     *
     * @param layout Grid layout used to iterate and locate edges.
     * @param signed_distance_at_nodes Node-centered signed-distance field.
     * @param cell_status Cell meshData used to detect ghost adjacency.
     * @param edge_status Edge meshData updated in place.
     */
    void tagEdges_(GridLayoutT const& layout,
                   node_double_field_type const& signed_distance_at_nodes,
                   cell_status_field_type const& cell_status,
                   edge_status_vecfield_type& edge_status) const
    {
        for (std::size_t dir = 0; dir < dim; ++dir)
        {
            layout.evalOnGhostBox(edge_status[dir], [&](auto... idx) {
                auto const local_edge = local_index_type{static_cast<std::uint32_t>(idx)...};
                auto const amr_edge   = fieldAMRIndex_(layout, edge_status[dir], local_edge);

                double phi_min = std::numeric_limits<double>::max();
                double phi_max = std::numeric_limits<double>::lowest();
                local_index_type node_idx{};
                forEachEdgeEndpoint_<0>(local_edge, dir, node_idx, [&](auto const& endpoint) {
                    auto const phi = signed_distance_at_nodes(endpoint);
                    phi_min        = std::min(phi_min, phi);
                    phi_max        = std::max(phi_max, phi);
                });

                if (isCut_(phi_min, phi_max, params_.cut_eps))
                {
                    edge_status[dir](local_edge) = toDouble(EdgeStatus::Cut);
                    return;
                }

                // Edges follow the same precedence as faces: cut first, then
                // ghost if they border the ghost-cell shell.
                if (edgeSurroundsGhostCell_(local_edge, dir, cell_status))
                {
                    edge_status[dir](local_edge) = toDouble(EdgeStatus::Ghost);
                    return;
                }

                auto const phi_edge = boundary_.signedDistance(
                    layout.fieldNodeCoordinates(edge_status[dir], amr_edge));
                edge_status[dir](local_edge)
                    = (phi_edge < -params_.inactive_eps) ? toDouble(EdgeStatus::Inactive) : toDouble(EdgeStatus::Fluid);
            });
        }
    }

    /**
     * @brief Validate that all fields use the centerings expected by the tagger.
     *
     * @param signed_distance_at_nodes Node-centered signed-distance field.
     * @param cell_status Cell-centered tag field.
     * @param face_status Face-centered tag vector field.
     * @param edge_status Edge-centered tag vector field.
     *
     * @throws std::runtime_error If one of the field centerings is inconsistent
     * with the tagger expectations.
     */
    static void validateCenterings_(node_double_field_type const& signed_distance_at_nodes,
                                    cell_status_field_type const& cell_status,
                                    face_status_vecfield_type const& face_status,
                                    edge_status_vecfield_type const& edge_status)
    {
        if (GridLayoutT::centering(signed_distance_at_nodes) != nodeCentering_())
            throw std::runtime_error("signed_distance_at_nodes has invalid centering");

        if (GridLayoutT::centering(cell_status) != cellCentering_())
            throw std::runtime_error("cell_status has invalid centering");

        auto const face_centering = GridLayoutT::centering(face_status);
        for (std::size_t dir = 0; dir < dim; ++dir)
            if (face_centering[dir] != faceCentering_(dir))
                throw std::runtime_error("face_status has invalid centering");

        auto const edge_centering = GridLayoutT::centering(edge_status);
        for (std::size_t dir = 0; dir < dim; ++dir)
            if (edge_centering[dir] != edgeCentering_(dir))
                throw std::runtime_error("edge_status has invalid centering");
    }
};

} // namespace PHARE::core

#endif
