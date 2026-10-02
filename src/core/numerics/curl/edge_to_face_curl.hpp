#ifndef PHARE_CORE_NUMERICS_CURL_EDGE_TO_FACE_CURL_HPP
#define PHARE_CORE_NUMERICS_CURL_EDGE_TO_FACE_CURL_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <array>
#include <cstddef>

namespace PHARE::core
{
// edge (E) centred scalar quantities, Ex/Ey/Ez share names in HybridQuantity and MHDQuantity
template<typename Scalar>
constexpr std::array<Scalar, 3> edge_quantities()
{
    return {Scalar::Ex, Scalar::Ey, Scalar::Ez};
}


/** @brief B = curl A, with A edge centred (E centring) and B face centred.
 *
 * Uses the same two-point staggered stencil (GridLayout::deriv) as Faraday, so the
 * discrete divergence of B vanishes to round-off. B is computed over its full ghost box:
 * each B component is dual in the directions it is differentiated in while its A source
 * is primal there, and primal/dual ghost widths are equal, so the A ghost box holds
 * every stencil node.
 */
template<typename GridLayout, typename AField, typename VecField>
void curl_edges_to_faces(GridLayout const& layout, std::array<AField, 3> const& A, VecField& B)
{
    constexpr auto dimension = GridLayout::dimension;
    static_assert(dimension > 1, "curl of a vector potential needs at least 2 dimensions");

    auto const& [Ax, Ay, Az] = A;
    auto& Bx                 = B(Component::X);
    auto& By                 = B(Component::Y);
    auto& Bz                 = B(Component::Z);

    layout.evalOnGhostBox(Bx, [&](auto const&... ijk) {
        if constexpr (dimension == 2)
            Bx(ijk...) = layout.template deriv<Direction::Y>(Az, {ijk...});
        else
            Bx(ijk...) = layout.template deriv<Direction::Y>(Az, {ijk...})
                         - layout.template deriv<Direction::Z>(Ay, {ijk...});
    });

    layout.evalOnGhostBox(By, [&](auto const&... ijk) {
        if constexpr (dimension == 2)
            By(ijk...) = -layout.template deriv<Direction::X>(Az, {ijk...});
        else
            By(ijk...) = layout.template deriv<Direction::Z>(Ax, {ijk...})
                         - layout.template deriv<Direction::X>(Az, {ijk...});
    });

    layout.evalOnGhostBox(Bz, [&](auto const&... ijk) {
        Bz(ijk...) = layout.template deriv<Direction::X>(Ay, {ijk...})
                     - layout.template deriv<Direction::Y>(Ax, {ijk...});
    });
}

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_CURL_EDGE_TO_FACE_CURL_HPP
