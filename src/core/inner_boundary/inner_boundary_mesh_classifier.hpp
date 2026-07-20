#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_CLASSIFIER_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_CLASSIFIER_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/inner_boundary/inner_boundary_geometry.hpp"
#include "core/inner_boundary/inner_boundary_mesh_data.hpp"
#include "core/numerics/interpolator/field_at_point.hpp"

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
 * @brief Classifies cells, faces, edges, and nodes around an embedded inner boundary.
 *
 * The classifier first evaluates the boundary signed distance on the node support,
 * classifies cells as fluid/cut/inactive, then grows a ghost-cell shell around
 * cut cells and finally propagates that geometry to all other element types.
 *
 * Element types are determined by their centering pattern.  For a dim-dimensional
 * simulation there are 2^dim distinct patterns:
 *  - all-dual  (n_primal == 0)         → cells
 *  - n_primal == 1                     → faces (one per direction, dim total)
 *  - 1 < n_primal < dim                → edges (one per direction, dim total; 3D only)
 *  - all-primal (n_primal == dim)      → nodes
 *
 * In 2D the four patterns are cell, face-X, face-Y, and node (e.g. Ez).
 * In 1D the two patterns are cell and node/face (which coincide in 1D).
 *
 * @tparam dim Spatial dimension.
 * @tparam GridLayoutT Grid layout type used to provide coordinates and iterate
 * local supports.
 * @tparam PhysicalQuantityT Quantity traits used to define node/cell/face/edge
 * field types.
 */
template<std::size_t dim, typename GridLayoutT, typename PhysicalQuantityT>
class InnerBoundaryMeshClassifier
{
public:
    using point_type              = InnerBoundaryGeometry<dim>::point_type;
    using local_index_type        = Point<std::uint32_t, dim>;
    using signed_local_index_type = Point<int, dim>;
    using mesh_data_type          = InnerBoundaryMeshData<dim, PhysicalQuantityT>;
    using field_type              = mesh_data_type::field_type;

    /**
     * @brief Runtime parameters controlling cut and ghost classification.
     *
     * @param nghosts Number of fluid-cell layers promoted to ghost around
     * cut cells.
     * @param cut_eps Tolerance used when deciding whether a support intersects
     * the boundary.
     * @param inactive_eps Tolerance used when deciding whether a non-cut support
     * lies inside the boundary.
     */
    struct Params
    {
        std::size_t nghosts = 1;
        double cut_eps      = 0.0;
        double inactive_eps = 0.0;

        /// True when the level under-resolves the boundary: the ghost shell is then grown only
        /// 1 cell deep and a degraded (first-order, ideal) flux is used in a band near the body.
        bool degraded = false;

        /// Full reconstruction reach of the solver (= the field ghost width). Used to size the
        /// degraded band even when @ref nghosts has been reduced to 1 for the shallow shell.
        std::size_t full_reach = 1;
    };

    /**
     * @brief Optional overrides for default parameter construction.
     *
     * Any unset entry keeps the mesh-derived default.
     */
    struct Overrides
    {
        std::optional<std::size_t> nghosts;
        std::optional<double> cut_eps;
        std::optional<double> inactive_eps;
    };

    /**
     * @brief Build a classifier with mesh-derived default tolerances.
     *
     * @param boundary Embedded boundary used for signed-distance queries.
     * @param layout Grid layout used to derive default tolerances from the mesh
     * spacing.
     * @param overrides Optional parameter overrides.
     * @return Configured classifier instance.
     */
    static InnerBoundaryMeshClassifier withDefaults(InnerBoundaryGeometry<dim> const& boundary,
                                                    GridLayoutT const& layout,
                                                    Overrides const& overrides = {})
    {
        auto const& dx    = layout.meshSize();
        auto const dx_min = *std::min_element(dx.begin(), dx.end());
        // Inner-boundary ghost shell / degraded-band reach is sized from the reconstruction
        // stencil (reconstruction_nghosts + 1) when the layout exposes it (MHD), NOT the field
        // ghost width (layout.nbrGhosts()). The field ghost width is inflated (even-rounded, with
        // extra layers for the J Laplacian and refinement formulas) and would over-grow the shell,
        // spuriously flipping well-resolved refined levels into the under-resolved/degraded path.
        // The reconstruction reach is what actually bounds how far the mirror-point stencil reads.
        // Layouts without a reconstruction stencil (e.g. Hybrid) fall back to the field ghosts.
        auto const nghosts = [&]() -> std::uint32_t {
            if constexpr (requires { GridLayoutT::implT::reconstruction_nghosts; })
                return GridLayoutT::implT::reconstruction_nghosts + 1;
            else
                return layout.nbrGhosts();
        }();

        // A level under-resolves the boundary when its characteristic length spans fewer than
        // (nghosts + 1) cells: then the ghost shell cannot be grown to the full stencil reach and
        // the mirror-point stencil would read into the body. On such levels the shell is grown
        // only 1 cell deep and a degraded (first-order, ideal) flux is used near the body, so the
        // stencil reaches only one cell into it. Resolved levels (e.g. finer levels, or a plane,
        // whose characteristicLength is infinite) are untouched.
        bool const underResolved
            = boundary.characteristicLength() < (static_cast<double>(nghosts) + 1.0) * dx_min;

        Params p;
        p.full_reach   = nghosts;
        p.degraded     = underResolved;
        p.nghosts      = overrides.nghosts.value_or(underResolved ? std::size_t{1} : nghosts);
        p.cut_eps      = overrides.cut_eps.value_or(1e-6 * dx_min);
        p.inactive_eps = overrides.inactive_eps.value_or(p.cut_eps);
        return InnerBoundaryMeshClassifier{boundary, p};
    }

    /**
     * @brief Construct a classifier with explicit parameters.
     *
     * @param boundary Embedded boundary used for signed-distance queries.
     * @param params Explicit cut/ghost classification parameters.
     */
    InnerBoundaryMeshClassifier(InnerBoundaryGeometry<dim> const& boundary, Params params)
        : boundary_{boundary}
        , params_{params}
    {
    }

    /**
     * @brief Fill node level-set values, classify all element types,
     *        and populate the precomputed ghost lists in @p meshData.
     *
     * @param layout Grid layout used to iterate and locate all supports.
     * @param meshData Bundle written by the classifier.
     */
    void operator()(GridLayoutT const& layout, mesh_data_type& meshData) const
    {
        validateCenterings_(meshData.signedDistanceAtNodes, meshData);
        fillNodePhi_(layout, meshData.signedDistanceAtNodes);

        auto& cell_status = meshData.cellStatusField();
        classifyCutInactiveAndFluidCells_(layout, meshData.signedDistanceAtNodes, cell_status);
        classifyGhostCells_(layout, cell_status);
        classifyNonCellElems_(layout, meshData.signedDistanceAtNodes, cell_status, meshData);
        populateGhostLists_(layout, meshData);
        populateDegradedLists_(layout, meshData);
    }

private:
    InnerBoundaryGeometry<dim> const& boundary_;
    Params params_;

    /**
     * @brief Decide whether a support is geometrically cut by the boundary.
     */
    static bool isCut_(double phi_min, double phi_max, double cut_eps)
    {
        if (phi_min < -cut_eps && phi_max > cut_eps)
            return true;
        return phi_min <= cut_eps && phi_max >= -cut_eps;
    }

    /**
     * @brief Check whether a signed local index lies inside a local array shape.
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
     */
    template<typename FieldT>
    static signed_local_index_type fieldAMRIndex_(GridLayoutT const& layout, FieldT const& field,
                                                  local_index_type const& local_idx)
    {
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

    /// All-primal centering (node-centered).
    static auto nodeCentering_()
    {
        std::array<QtyCentering, dim> centering{};
        centering.fill(QtyCentering::primal);
        return centering;
    }

    /// All-dual centering (cell-centered).
    static auto cellCentering_()
    {
        std::array<QtyCentering, dim> centering{};
        centering.fill(QtyCentering::dual);
        return centering;
    }

    // -------------------------------------------------------------------------
    //                      support-node and adjacent-cell iteration
    // -------------------------------------------------------------------------

    /**
     * @brief Visit all nodes in the support of a mesh element.
     *
     * The centering of the element drives the traversal: a **dual** direction
     * contributes both of the two bounding nodes (indices @c elem[i] and
     * @c elem[i]+1), while a **primal** direction contributes a single node
     * (index @c elem[i]).  The total number of visits is 2^(number of dual
     * directions).
     *
     * Examples for dim=2:
     *  - all-dual (cell): visits all 4 corner nodes.
     *  - one-primal/one-dual (face): visits the 2 nodes bounding the face.
     *  - all-primal (node, e.g. Ez in MHD): visits the element itself.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachSupportNode_(local_index_type const& elem,
                                    std::array<QtyCentering, dim> const& centering,
                                    local_index_type& node, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(node);
        }
        else if (centering[iDim] == QtyCentering::dual)
        {
            node[iDim] = elem[iDim];
            forEachSupportNode_<iDim + 1>(elem, centering, node, std::forward<Fn>(fn));
            node[iDim] = elem[iDim] + 1;
            forEachSupportNode_<iDim + 1>(elem, centering, node, std::forward<Fn>(fn));
        }
        else // primal
        {
            node[iDim] = elem[iDim];
            forEachSupportNode_<iDim + 1>(elem, centering, node, std::forward<Fn>(fn));
        }
    }

    /**
     * @brief Visit all cells that contain (are adjacent to) a mesh element.
     *
     * The rule is the complement of @c forEachSupportNode_: a **primal**
     * direction contributes two adjacent cells (at @c elem[i]-1 and @c elem[i]),
     * while a **dual** direction contributes a single cell (at @c elem[i]).
     * The total number of visits is 2^(number of primal directions).
     *
     * Examples for dim=2:
     *  - all-dual (cell): visits the cell itself (trivially contained in itself).
     *  - one-primal/one-dual (face): visits the 2 cells on either side of the face.
     *  - all-primal (node, e.g. Ez): visits all 4 surrounding cells.
     */
    template<std::size_t iDim, typename Fn>
    static void forEachAdjacentCell_(local_index_type const& elem,
                                     std::array<QtyCentering, dim> const& centering,
                                     signed_local_index_type& cell, Fn&& fn)
    {
        if constexpr (iDim == dim)
        {
            fn(cell);
        }
        else if (centering[iDim] == QtyCentering::primal)
        {
            cell[iDim] = static_cast<int>(elem[iDim]) - 1;
            forEachAdjacentCell_<iDim + 1>(elem, centering, cell, std::forward<Fn>(fn));
            cell[iDim] = static_cast<int>(elem[iDim]);
            forEachAdjacentCell_<iDim + 1>(elem, centering, cell, std::forward<Fn>(fn));
        }
        else // dual
        {
            cell[iDim] = static_cast<int>(elem[iDim]);
            forEachAdjacentCell_<iDim + 1>(elem, centering, cell, std::forward<Fn>(fn));
        }
    }

    // -------------------------------------------------------------------------
    //                          ghost-adjacency helper
    // -------------------------------------------------------------------------

    /**
     * @brief Return true iff any cell adjacent to @p elem (given its @p centering)
     * has Fluid or Cut status.
     */
    bool elemSurroundsFluidOrCutCell_(local_index_type const& elem,
                                      std::array<QtyCentering, dim> const& centering,
                                      field_type const& cell_status) const
    {
        bool found            = false;
        auto const cell_shape = cell_status.shape();
        signed_local_index_type cell{};
        forEachAdjacentCell_<0>(elem, centering, cell, [&](auto const& adjacent_cell) {
            if (found || !inBounds_(adjacent_cell, cell_shape))
                return;
            auto const s = cell_status(asLocal_(adjacent_cell));
            if (s == toDouble(ElemStatus::Fluid) || s == toDouble(ElemStatus::Cut))
                found = true;
        });
        return found;
    }

    /**
     * @brief Return true iff any cell adjacent to @p elem (given its @p centering)
     * has Ghost status.
     *
     * Delegates to forEachAdjacentCell_, which visits 2^(n_primal) cells based on
     * the centering — no special-casing for faces, edges, or nodes required.
     */
    bool elemSurroundsGhostCell_(local_index_type const& elem,
                                 std::array<QtyCentering, dim> const& centering,
                                 field_type const& cell_status) const
    {
        bool has_ghost        = false;
        auto const cell_shape = cell_status.shape();
        signed_local_index_type cell{};
        forEachAdjacentCell_<0>(elem, centering, cell, [&](auto const& adjacent_cell) {
            if (has_ghost || !inBounds_(adjacent_cell, cell_shape))
                return;
            if (cell_status(asLocal_(adjacent_cell)) == toDouble(ElemStatus::Ghost))
                has_ghost = true;
        });
        return has_ghost;
    }

    // -------------------------------------------------------------------------
    //                           classification passes
    // -------------------------------------------------------------------------

    void fillNodePhi_(GridLayoutT const& layout, field_type& signed_distance_at_nodes) const
    {
        layout.evalOnGhostBox(signed_distance_at_nodes, [&](auto... idx) {
            auto const local_point = local_index_type{static_cast<std::uint32_t>(idx)...};
            auto const amr_point   = fieldAMRIndex_(layout, signed_distance_at_nodes, local_point);
            signed_distance_at_nodes(local_point) = boundary_.signedDistance(
                layout.fieldNodeCoordinates(signed_distance_at_nodes, amr_point));
        });
    }

    void classifyGhostCells_(GridLayoutT const& layout, field_type& cell_status) const
    {
        if (params_.nghosts == 0)
            return;

        auto const cell_shape = cell_status.shape();
        std::vector<signed_local_index_type> frontier;

        layout.evalOnGhostBox(cell_status, [&](auto... idx) {
            auto const cell = local_index_type{static_cast<std::uint32_t>(idx)...};
            if (cell_status(cell) == toDouble(ElemStatus::Cut))
                frontier.push_back(asSigned_(cell));
        });

        for (std::size_t layer = 0; layer < params_.nghosts && !frontier.empty(); ++layer)
        {
            std::vector<signed_local_index_type> next_frontier;
            next_frontier.reserve(frontier.size() * 2 * dim);

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

                        if (cell_status(local_neigh) != toDouble(ElemStatus::Inactive))
                            continue;

                        cell_status(local_neigh) = toDouble(ElemStatus::Ghost);
                        next_frontier.push_back(neigh);
                    }
                }
            }

            frontier = std::move(next_frontier);
        }
    }

    void classifyCutInactiveAndFluidCells_(GridLayoutT const& layout,
                                           field_type const& signed_distance_at_nodes,
                                           field_type& cell_status) const
    {
        auto const centering = cellCentering_();
        layout.evalOnGhostBox(cell_status, [&](auto... idx) {
            auto const local_cell = local_index_type{static_cast<std::uint32_t>(idx)...};
            auto const amr_cell   = fieldAMRIndex_(layout, cell_status, local_cell);

            double phi_min = std::numeric_limits<double>::max();
            double phi_max = std::numeric_limits<double>::lowest();
            local_index_type node_idx{};
            forEachSupportNode_<0>(local_cell, centering, node_idx, [&](auto const& node) {
                auto const phi = signed_distance_at_nodes(node);
                phi_min        = std::min(phi_min, phi);
                phi_max        = std::max(phi_max, phi);
            });

            if (isCut_(phi_min, phi_max, params_.cut_eps))
            {
                cell_status(local_cell) = toDouble(ElemStatus::Cut);
                return;
            }

            auto const phi_center
                = boundary_.signedDistance(layout.cellCenteredCoordinates(amr_cell));
            cell_status(local_cell) = (phi_center < -params_.inactive_eps)
                                          ? toDouble(ElemStatus::Inactive)
                                          : toDouble(ElemStatus::Fluid);
        });
    }

    /**
     * @brief Classify all non-cell element types (faces, edges, nodes).
     *
     * Iterates over all 2^dim centering patterns and skips the all-dual (cell)
     * pattern which is handled by a dedicated pass. For every other pattern the
     * same logic applies regardless of whether the element is a face, an edge,
     * or a node: forEachSupportNode_ and forEachAdjacentCell_ adapt automatically
     * to the centering.
     */
    void classifyNonCellElems_(GridLayoutT const& layout,
                               field_type const& signed_distance_at_nodes,
                               field_type const& cell_status, mesh_data_type& meshData) const
    {
        for (std::size_t idx = 0; idx < mesh_data_type::num_elem_types; ++idx)
        {
            auto const centering = mesh_data_type::idxToCentering(idx);
            auto const n_dual    = static_cast<std::size_t>(
                std::count(centering.begin(), centering.end(), QtyCentering::dual));

            if (n_dual == dim)
                continue; // cell — handled by classifyCutInactiveAndFluidCells_

            auto& elem_field = meshData.elemStatus[idx];
            layout.evalOnGhostBox(elem_field, [&](auto... local_idx_args) {
                auto const local_elem
                    = local_index_type{static_cast<std::uint32_t>(local_idx_args)...};
                auto const amr_elem = fieldAMRIndex_(layout, elem_field, local_elem);

                double phi_min = std::numeric_limits<double>::max();
                double phi_max = std::numeric_limits<double>::lowest();
                local_index_type node_idx{};
                forEachSupportNode_<0>(local_elem, centering, node_idx,
                                       [&](auto const& support_node) {
                                           auto const phi = signed_distance_at_nodes(support_node);
                                           phi_min        = std::min(phi_min, phi);
                                           phi_max        = std::max(phi_max, phi);
                                       });

                if (isCut_(phi_min, phi_max, params_.cut_eps))
                {
                    elem_field(local_elem) = toDouble(ElemStatus::Cut);
                    return;
                }

                // let too close elements from the boundary be active to avoid complicated
                // interpolation at the symmetric location
                if (elemSurroundsFluidOrCutCell_(local_elem, centering, cell_status))
                {
                    elem_field(local_elem) = toDouble(ElemStatus::Fluid);
                    return;
                }

                if (elemSurroundsGhostCell_(local_elem, centering, cell_status))
                {
                    elem_field(local_elem) = toDouble(ElemStatus::Ghost);
                    return;
                }

                auto const phi_elem
                    = boundary_.signedDistance(layout.fieldNodeCoordinates(elem_field, amr_elem));
                elem_field(local_elem) = (phi_elem < -params_.inactive_eps)
                                             ? toDouble(ElemStatus::Inactive)
                                             : toDouble(ElemStatus::Fluid);
            });
        }
    }

    /**
     * @brief Scan all status fields and populate precomputed ghost lists in @p meshData.
     *
     * Iterates over all 2^dim centering patterns. For each ghost element the
     * precomputed data are:
     *  - its local array index,
     *  - the physical position of its symmetric (mirror) point in the fluid,
     *  - the outward boundary normal at the element centre.
     */
    void populateGhostLists_(GridLayoutT const& layout, mesh_data_type& meshData) const
    {
        for (auto& vec : meshData.ghostElemsData)
            vec.clear();

        for (std::size_t idx = 0; idx < mesh_data_type::num_elem_types; ++idx)
        {
            auto const& status_field = meshData.elemStatus[idx];
            auto& ghost_list         = meshData.ghostElemsData[idx];
            auto const centerings    = mesh_data_type::idxToCentering(idx);

            layout.evalOnGhostBox(status_field, [&](auto... local_idx_args) {
                auto const local = local_index_type{static_cast<std::uint32_t>(local_idx_args)...};
                if (status_field(local) != toDouble(ElemStatus::Ghost))
                    return;

                auto const amr    = fieldAMRIndex_(layout, status_field, local);
                auto const pos    = layout.fieldNodeCoordinates(status_field, amr);
                auto const mirror = boundary_.symmetric(pos);
                ghost_list.push_back(
                    {local, mirror, boundary_.normal(pos),
                     FieldAtPoint<dim, 1>::pointIsInterpolable(layout, mirror, centerings)});
            });
        }
    }

    /**
     * @brief Populate the per-centering degraded-element lists.
     *
     * On under-resolved levels, flags every mesh element whose full reconstruction stencil
     * (reaching @c full_reach cells) could touch the body: the flux and CT correction passes
     * recompute those elements with a first-order, ideal scheme. The flag set is built by
     * growing a @c full_reach -cell ball outward from the Cut/Ghost cells over the cell grid,
     * then marking any element adjacent to a ball cell. On resolved levels the lists stay empty.
     */
    void populateDegradedLists_(GridLayoutT const& layout, mesh_data_type& meshData) const
    {
        for (auto& vec : meshData.degradedElemsData)
            vec.clear();

        if (!params_.degraded)
            return;

        auto const& cell_status = meshData.cellStatusField();
        auto const cell_shape   = cell_status.shape();

        std::size_t total = 1;
        for (std::size_t d = 0; d < dim; ++d)
            total *= cell_shape[d];

        auto flatten = [&](signed_local_index_type const& c) {
            std::size_t flat = 0;
            for (std::size_t d = 0; d < dim; ++d)
                flat = flat * cell_shape[d] + static_cast<std::size_t>(c[d]);
            return flat;
        };

        // Seed the ball with all Cut/Ghost cells, then grow full_reach layers outward into any
        // in-bounds cell (Fluid as well as Inactive — the stencil reach is geometric).
        std::vector<char> nearBody(total, 0);
        std::vector<signed_local_index_type> frontier;
        layout.evalOnGhostBox(cell_status, [&](auto... idx) {
            auto const cell = local_index_type{static_cast<std::uint32_t>(idx)...};
            auto const s    = cell_status(cell);
            if (s == toDouble(ElemStatus::Cut) || s == toDouble(ElemStatus::Ghost))
            {
                auto const signed_cell  = asSigned_(cell);
                nearBody[flatten(signed_cell)] = 1;
                frontier.push_back(signed_cell);
            }
        });

        for (std::size_t layer = 0; layer < params_.full_reach && !frontier.empty(); ++layer)
        {
            std::vector<signed_local_index_type> next_frontier;
            for (auto const& source : frontier)
            {
                for (std::size_t d = 0; d < dim; ++d)
                {
                    for (int sgn : {-1, 1})
                    {
                        auto neigh = source;
                        neigh[d] += sgn;
                        if (!inBounds_(neigh, cell_shape) || nearBody[flatten(neigh)])
                            continue;
                        nearBody[flatten(neigh)] = 1;
                        next_frontier.push_back(neigh);
                    }
                }
            }
            frontier = std::move(next_frontier);
        }

        // Mark every element (any centering) that touches a near-body cell.
        // Iterate the *physical* box only (not the ghost box): the flux/E degrade passes
        // reconstruct at each listed element with a first-order stencil reading the
        // cell-centered state at idx and previous(idx). That stencil is in-bounds exactly on
        // the physical box (the same domain the normal flux/E passes use via evalOnBox); a
        // ghost-region face/edge would read the cell-centered field out of its allocated
        // array (face fields have one more point than cell fields along a primal axis),
        // causing a segfault. Ghost-region elements are never consumed by the solver anyway.
        for (std::size_t idx = 0; idx < mesh_data_type::num_elem_types; ++idx)
        {
            auto const centering     = mesh_data_type::idxToCentering(idx);
            auto const& status_field = meshData.elemStatus[idx];
            auto& degraded_list      = meshData.degradedElemsData[idx];

            layout.evalOnBox(status_field, [&](auto... local_idx_args) {
                auto const local_elem
                    = local_index_type{static_cast<std::uint32_t>(local_idx_args)...};

                bool degraded = false;
                signed_local_index_type cell{};
                forEachAdjacentCell_<0>(local_elem, centering, cell,
                                        [&](auto const& adjacent_cell) {
                                            if (degraded || !inBounds_(adjacent_cell, cell_shape))
                                                return;
                                            if (nearBody[flatten(adjacent_cell)])
                                                degraded = true;
                                        });

                if (degraded)
                    degraded_list.push_back({local_elem, {}, {}, false});
            });
        }
    }

    /**
     * @brief Validate that all fields carry the centerings expected by the classifier.
     *
     * @throws std::runtime_error on any centering mismatch.
     */
    static void validateCenterings_(field_type const& signed_distance_at_nodes,
                                    mesh_data_type const& meshData)
    {
        if (GridLayoutT::centering(signed_distance_at_nodes) != nodeCentering_())
            throw std::runtime_error("signed_distance_at_nodes has invalid centering");

        for (std::size_t idx = 0; idx < mesh_data_type::num_elem_types; ++idx)
        {
            auto const expected = mesh_data_type::idxToCentering(idx);
            if (GridLayoutT::centering(meshData.elemStatus[idx]) != expected)
                throw std::runtime_error("elemStatus entry at index " + std::to_string(idx)
                                         + " has invalid centering");
        }
    }
};

} // namespace PHARE::core

#endif
