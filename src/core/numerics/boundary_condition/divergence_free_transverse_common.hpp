#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_COMMON_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_COMMON_HPP

#include "core/data/grid/gridlayoutdefs.hpp"

#include <cstddef>
#include <tuple>
#include <type_traits>

namespace PHARE::core
{
/**
 * @brief Enforce the discrete divergence-free condition on the boundary-normal magnetic component,
 * given the transverse components already filled in the ghost layer.
 *
 * Iterates the ghost cells closest-to-farthest from the physical domain (forward for an upper
 * boundary, reverse for a lower one) so each normal update reads an already-set donor.
 *
 * Shared by the transverse-Neumann and transverse-Dirichlet divergence-free conditions.
 *
 * @param fields         Tuple of the three (co-centred) magnetic-field component fields.
 * @param iNormal        Index of the boundary-normal component.
 * @param side           Which side of the domain the boundary is on.
 * @param gridLayout     Grid layout (mesh spacings).
 * @param localGhostBox  Local ghost box swept for the normal-component update.
 */
template<std::size_t dimension, typename FieldTuple, typename GridLayoutT, typename BoxT>
void applyDivergenceFreeNormalComponent(FieldTuple& fields, std::size_t const iNormal,
                                        Side const side, GridLayoutT const& gridLayout,
                                        BoxT const& localGhostBox)
{
    using _index     = Point<std::uint32_t, dimension>;
    using field_type = std::tuple_element_t<0, std::remove_reference_t<FieldTuple>>;

    field_type& nField = [&]() -> field_type& {
        switch (iNormal)
        {
            case 0: return std::get<0>(fields);
            case 1: return std::get<1>(fields);
            default: return std::get<2>(fields);
        }
    }();

    auto apply_loop = [&](auto begin, auto end) {
        for (auto it = begin; it != end; ++it)
        {
            _index const& index = *it;

            double transverseDiv = 0.0;
            for_N<dimension>([&](auto iTransverse) {
                if (static_cast<std::size_t>(iTransverse) != iNormal)
                {
                    field_type& tField       = std::get<iTransverse>(fields);
                    _index const upper_index = index.shift(iTransverse, 1);
                    double const invDxT      = gridLayout.inverseMeshSize(
                        static_cast<Direction>(static_cast<std::uint32_t>(iTransverse)));
                    transverseDiv += (tField(upper_index) - tField(index)) * invDxT;
                }
            });

            double const dxN = gridLayout.meshSize()[iNormal];
            if (side == Side::Upper)
            {
                _index const index_to_set      = index.shift(iNormal, 1);
                _index const index_already_set = index;
                nField(index_to_set)           = nField(index_already_set) - dxN * transverseDiv;
            }
            else
            {
                _index const index_to_set      = index;
                _index const index_already_set = index.shift(iNormal, 1);
                nField(index_to_set)           = nField(index_already_set) + dxN * transverseDiv;
            }
        }
    };

    if (side == Side::Upper)
        apply_loop(localGhostBox.begin(), localGhostBox.end());
    else
        apply_loop(localGhostBox.rbegin(), localGhostBox.rend());
}

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_COMMON_HPP
