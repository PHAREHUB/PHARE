#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_DATA_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_DATA_HPP

#include "core/def.hpp"
#include "core/data/field/field.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/inner_boundary/ghost_elem_pack.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"
#include "core/utilities/point/point.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace PHARE::core
{
/**
 * @brief Bundle of node level-set values and mesh status data around an inner boundary.
 *
 * Status fields and ghost data are indexed by a binary encoding of the element centering:
 * bit i of the index is 1 when direction i is dual, 0 when primal. This gives 2^dim
 * distinct element types stored in fixed-size arrays.
 *
 * In a 2D Yee grid the four types are, by index:
 *   0 = (P,P) node-centered    (e.g. Ez)
 *   1 = (D,P) face-Y-centered  (e.g. Bx)
 *   2 = (P,D) face-X-centered  (e.g. By)
 *   3 = (D,D) cell-centered    (e.g. rho, Bz)
 *
 * Use @ref centeringToIdx to convert a centering array to its index and
 * @ref getStatusFieldFromCentering / @ref getGhostDataFromCentering for typed access.
 *
 * @tparam dim Spatial dimension.
 * @tparam PhysicalQuantityT Quantity traits providing field scalar types and centering info.
 */
template<std::size_t dim, typename PhysicalQuantityT>
struct InnerBoundaryMeshData
{
    using field_type           = Field<dim, typename PhysicalQuantityT::Scalar, double>;
    using vecfield_type        = VecField<field_type, PhysicalQuantityT>;
    using ghost_elem_data_type = GhostElemData<dim>;

    /// Number of distinct element types: one per centering pattern (2^dim).
    static constexpr std::size_t num_elem_types = (std::size_t{1} << dim);


    // -------------------------------------------------------------------------
    //                        centering ↔ index helpers
    // -------------------------------------------------------------------------

    /**
     * @brief Encode a centering pattern as an array index.
     *
     * Bit i of the returned index is 1 when direction i is @c dual.
     * This gives a unique index in [0, 2^dim) for every centering pattern.
     *
     * Examples (dim=2):  (P,P)→0,  (D,P)→1,  (P,D)→2,  (D,D)→3
     */
    static constexpr std::size_t
    centeringToIdx(std::array<QtyCentering, dim> const& c) noexcept
    {
        std::size_t idx = 0;
        for (std::size_t i = 0; i < dim; ++i)
            if (c[i] == QtyCentering::dual)
                idx |= (std::size_t{1} << i);
        return idx;
    }

    /**
     * @brief Decode an array index back to a centering pattern (inverse of centeringToIdx).
     */
    static constexpr auto idxToCentering(std::size_t idx) noexcept
    {
        std::array<QtyCentering, dim> c{};
        for (std::size_t i = 0; i < dim; ++i)
            c[i] = ((idx >> i) & 1u) ? QtyCentering::dual : QtyCentering::primal;
        return c;
    }

    /**
     * @brief Return the PhysicalQuantityT::Scalar whose centering matches @p centering.
     *
     * Mapping rule:
     *   - all dual   → CellCentered
     *   - all primal → NodeCentered
     *   - 1 primal   → FaceCentered[primal_dir]
     *   - otherwise  → EdgeCentered[dual_dir]  (exactly one dual direction in this case)
     */
    static auto scalarFromCentering(std::array<QtyCentering, dim> const& centering)
    {
        auto const n_primal = static_cast<std::size_t>(
            std::count(centering.begin(), centering.end(), QtyCentering::primal));

        if (n_primal == 0)
            return PhysicalQuantityT::Scalar::CellCentered;

        if (n_primal == dim)
            return PhysicalQuantityT::Scalar::NodeCentered;

        if (n_primal == 1)
        {
            auto const dir = static_cast<std::size_t>(std::distance(
                centering.begin(),
                std::find(centering.begin(), centering.end(), QtyCentering::primal)));
            return PhysicalQuantityT::FaceCentered()[dir];
        }

        // 1 < n_primal < dim: exactly one dual direction
        auto const dir = static_cast<std::size_t>(std::distance(
            centering.begin(),
            std::find(centering.begin(), centering.end(), QtyCentering::dual)));
        return PhysicalQuantityT::EdgeCentered()[dir];
    }


    // -------------------------------------------------------------------------
    //                              data members
    // -------------------------------------------------------------------------

    field_type signedDistanceAtNodes; ///< Signed distance to the boundary at nodes.

    /**
     * @brief Per-centering status fields, indexed by @ref centeringToIdx.
     *
     * Contains one status field per distinct centering pattern (2^dim entries).
     */
    std::array<field_type, num_elem_types> elemStatus;

    /**
     * @brief Per-centering ghost element lists, indexed by @ref centeringToIdx.
     *
     * Populated by InnerBoundaryMeshClassifier, iterated by the BC applier every time step.
     *
     * Backed by a per-patch GhostElemPatchData allocated by SAMRAI; rebound on each
     * setOnPatch via ResourcesManager.
     */
    GhostElemPack<dim> ghostElemsData;

    /**
     * @brief Per-centering lists of elements needing a degraded (first-order, ideal) flux/E
     * recompute, indexed by @ref centeringToIdx.
     *
     * Populated by the classifier only on levels that under-resolve the boundary
     * (characteristicLength/dx below the ghost-shell requirement): an element is listed when its
     * full reconstruction stencil would reach a Cut/Ghost cell. Empty on resolved levels.
     *
     * Reuses the GhostElemPack backing — only the @c index field of each entry is meaningful here
     * (mirror/normal are unused). The flux and CT correction passes iterate these lists to
     * overwrite the high-order result with a first-order, ideal computation.
     */
    GhostElemPack<dim> degradedElemsData;


    // -------------------------------------------------------------------------
    //                              constructor
    // -------------------------------------------------------------------------

    explicit InnerBoundaryMeshData(std::string boundaryName = "")
        : signedDistanceAtNodes{boundaryName + "_signed_distance",
                                PhysicalQuantityT::Scalar::NodeCentered}
        , elemStatus{makeElemStatus_(boundaryName, std::make_index_sequence<num_elem_types>{})}
        , ghostElemsData{boundaryName + "_ghost_elems"}
        , degradedElemsData{boundaryName + "_degraded_elems"}
    {
    }


    // -------------------------------------------------------------------------
    //                              accessors
    // -------------------------------------------------------------------------

    /**
     * @brief Return the status field for the mesh element type matching @p centering.
     */
    field_type& getStatusFieldFromCentering(std::array<QtyCentering, dim> const& centering)
    {
        return elemStatus[centeringToIdx(centering)];
    }

    field_type const&
    getStatusFieldFromCentering(std::array<QtyCentering, dim> const& centering) const
    {
        return elemStatus[centeringToIdx(centering)];
    }

    /**
     * @brief Return the cell-centered status field (all-dual centering).
     *
     * Convenience accessor — equivalent to getStatusFieldFromCentering with
     * all-dual centering, which is always elemStatus.back().
     */
    field_type& cellStatusField() { return elemStatus.back(); }
    field_type const& cellStatusField() const { return elemStatus.back(); }

    /**
     * @brief Return a tuple of status field references for a VecField component triplet.
     *
     * @param centerings Three per-dimension centering arrays, one per component.
     */
    auto getStatusFieldsFromCenterings(
        std::array<std::array<QtyCentering, dim>, 3> const& centerings)
        -> std::tuple<field_type&, field_type&, field_type&>
    {
        return std::forward_as_tuple(getStatusFieldFromCentering(centerings[0]),
                                     getStatusFieldFromCentering(centerings[1]),
                                     getStatusFieldFromCentering(centerings[2]));
    }

    auto getStatusFieldsFromCenterings(
        std::array<std::array<QtyCentering, dim>, 3> const& centerings) const
        -> std::tuple<field_type const&, field_type const&, field_type const&>
    {
        return std::forward_as_tuple(getStatusFieldFromCentering(centerings[0]),
                                     getStatusFieldFromCentering(centerings[1]),
                                     getStatusFieldFromCentering(centerings[2]));
    }

    /**
     * @brief Return the ghost data list for the mesh element type matching @p centering.
     */
    std::vector<ghost_elem_data_type>&
    getGhostDataFromCentering(std::array<QtyCentering, dim> const& centering)
    {
        return ghostElemsData[centeringToIdx(centering)];
    }

    std::vector<ghost_elem_data_type> const&
    getGhostDataFromCentering(std::array<QtyCentering, dim> const& centering) const
    {
        return ghostElemsData[centeringToIdx(centering)];
    }

    /**
     * @brief Return the degraded-element list for the mesh element type matching @p centering.
     */
    std::vector<ghost_elem_data_type>&
    getDegradedElemsFromCentering(std::array<QtyCentering, dim> const& centering)
    {
        return degradedElemsData[centeringToIdx(centering)];
    }

    std::vector<ghost_elem_data_type> const&
    getDegradedElemsFromCentering(std::array<QtyCentering, dim> const& centering) const
    {
        return degradedElemsData[centeringToIdx(centering)];
    }

    /**
     * @brief Return a tuple of ghost data list references for a VecField component triplet.
     *
     * @param centerings Three per-dimension centering arrays, one per component.
     */
    auto getGhostDataFromCenterings(
        std::array<std::array<QtyCentering, dim>, 3> const& centerings)
        -> std::tuple<std::vector<ghost_elem_data_type>&, std::vector<ghost_elem_data_type>&,
                      std::vector<ghost_elem_data_type>&>
    {
        return std::forward_as_tuple(getGhostDataFromCentering(centerings[0]),
                                     getGhostDataFromCentering(centerings[1]),
                                     getGhostDataFromCentering(centerings[2]));
    }

    auto getGhostDataFromCenterings(
        std::array<std::array<QtyCentering, dim>, 3> const& centerings) const
        -> std::tuple<std::vector<ghost_elem_data_type> const&,
                      std::vector<ghost_elem_data_type> const&,
                      std::vector<ghost_elem_data_type> const&>
    {
        return std::forward_as_tuple(getGhostDataFromCentering(centerings[0]),
                                     getGhostDataFromCentering(centerings[1]),
                                     getGhostDataFromCentering(centerings[2]));
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return std::apply(
            [this](auto const&... fields) {
                return isUsable(signedDistanceAtNodes, fields..., ghostElemsData, degradedElemsData);
            },
            elemStatus);
    }

    NO_DISCARD bool isSettable() const
    {
        return std::apply(
            [this](auto const&... fields) {
                return isSettable(signedDistanceAtNodes, fields..., ghostElemsData,
                                  degradedElemsData);
            },
            elemStatus);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::apply(
            [this](auto&... fields) {
                return std::forward_as_tuple(signedDistanceAtNodes, fields..., ghostElemsData,
                                             degradedElemsData);
            },
            elemStatus);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::apply(
            [this](auto const&... fields) {
                return std::forward_as_tuple(signedDistanceAtNodes, fields..., ghostElemsData,
                                             degradedElemsData);
            },
            elemStatus);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------


private:
    static std::string centeringToString_(std::array<QtyCentering, dim> const& c)
    {
        std::string s;
        for (auto const& qc : c)
            s += (qc == QtyCentering::primal ? 'P' : 'D');
        return s;
    }

    template<std::size_t... Is>
    static auto makeElemStatus_(std::string const& name, std::index_sequence<Is...>)
    {
        return std::array<field_type, sizeof...(Is)>{
            [&](std::size_t i) -> field_type {
                auto const c = idxToCentering(i);
                return field_type{name + "_status_" + centeringToString_(c),
                                  scalarFromCentering(c)};
            }(Is)...
        };
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_MESH_DATA_HPP
