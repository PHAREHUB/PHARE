#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_NEUMANN_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_NEUMANN_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/divergence_free_transverse_common.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include <algorithm>
#include <cstddef>

namespace PHARE::core
{
/**
 * @brief Boundary condition for the magnetic field B that enforces zero normal derivative on the
 * transverse components and sets the normal component so that the numerical divergence of B is
 * zero.
 *
 * On the transverse components the ghost value mirrors the first interior value (zero normal
 * gradient): B(index) = B(mirror). The normal component is then set so that div B = 0.
 *
 * @warning Only valid for vector fields with the same centering as the magnetic field.
 *
 * @tparam VecFieldT Type of the vector field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename VecFieldT, typename GridLayoutT>
class FieldBoundaryConditionDivergenceFreeTransverseNeumann
    : public IFieldBoundaryCondition<VecFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<VecFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static_assert(
        N == 3,
        "Divergence-free transverse Neumann boundary condition only applies to vector fields.");

    FieldBoundaryConditionDivergenceFreeTransverseNeumann() = default;

    FieldBoundaryConditionDivergenceFreeTransverseNeumann(
        FieldBoundaryConditionDivergenceFreeTransverseNeumann const&) = default;
    FieldBoundaryConditionDivergenceFreeTransverseNeumann&
    operator=(FieldBoundaryConditionDivergenceFreeTransverseNeumann const&) = default;
    FieldBoundaryConditionDivergenceFreeTransverseNeumann(
        FieldBoundaryConditionDivergenceFreeTransverseNeumann&&) = default;
    FieldBoundaryConditionDivergenceFreeTransverseNeumann&
    operator=(FieldBoundaryConditionDivergenceFreeTransverseNeumann&&) = default;

    virtual ~FieldBoundaryConditionDivergenceFreeTransverseNeumann() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::DivergenceFreeTransverseNeumann;
    }

    void apply(VecFieldT& vecField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               [[maybe_unused]] double const time) override
    {
        Direction const direction = getDirection(boundaryLocation);
        Side const side           = getSide(boundaryLocation);
        std::size_t const iNormal = static_cast<std::size_t>(direction);

        auto fields = vecField.components();

        assert(gridLayout.centering(vecField) == gridLayout.centering(tensor_quantity_type::B));

        // transverse components: zero normal gradient, B(index) = B(mirror)
        for_N<N>([&](auto iTransverse) {
            if (static_cast<std::size_t>(iTransverse) != iNormal)
            {
                field_type& Bc = std::get<iTransverse>(fields);

                // iNormal < dimension always holds: the std::min is a no-op that only lets GCC
                // prove the index is in bounds, otherwise the 1D instantiation (centering array of
                // size 1) triggers a false-positive -Werror=array-bounds in optimized builds.
                QtyCentering const centering = GridLayoutT::centering(
                    Bc.physicalQuantity())[std::min(iNormal, dimension - 1)];
                auto fieldBox = gridLayout.toFieldBox(localGhostBox, Bc.physicalQuantity());

                for (auto const& index : fieldBox)
                {
                    auto const mirrorIndex
                        = gridLayout.boundaryMirrored(direction, side, centering, index);
                    Bc(index) = Bc(mirrorIndex);
                }
            }
        });

        // set the normal component so the discrete divergence of B is zero, given the transverse
        // ghosts filled above (shared with the transverse-Dirichlet condition).
        applyDivergenceFreeNormalComponent<dimension>(fields, iNormal, side, gridLayout,
                                                      localGhostBox);
    }

}; // class FieldBoundaryConditionDivergenceFreeTransverseNeumann

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_NEUMANN_HPP
