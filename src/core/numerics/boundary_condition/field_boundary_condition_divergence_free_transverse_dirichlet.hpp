#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_DIRICHLET_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_DIRICHLET_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/divergence_free_transverse_common.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_dirichlet.hpp"

#include <array>
#include <cstddef>

namespace PHARE::core
{
/**
 * @brief Boundary condition for vector fields that imposes a value on tangential
 * components and sets the normal component so that numerical divergence is zero.
 *
 * Each tangential component is delegated to a scalar @c FieldBoundaryConditionDirichlet;
 * the normal component is then recomputed from the tangential ghost values so that the
 * discrete divergence stays zero.
 *
 * @warning Only valid for vector fields with the same centering as the magnetic field.
 *
 * @tparam VecFieldT Type of the vector field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename VecFieldT, typename GridLayoutT>
class FieldBoundaryConditionDivergenceFreeTransverseDirichlet
    : public IFieldBoundaryCondition<VecFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<VecFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;
    using value_type           = field_type::value_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static_assert(N == 3,
                  "Divergence-free transverse Dirichlet boundary condition only applies to vector "
                  "fields.");

    FieldBoundaryConditionDivergenceFreeTransverseDirichlet() = default;

    FieldBoundaryConditionDivergenceFreeTransverseDirichlet(value_type value,
                                                            DirichletExtrapolation extrapolation)
    {
        for (std::size_t i = 0; i < N; ++i)
            scalar_dirichlet_conditions_[i] = _scalar_dirichlet_bc_type{value, extrapolation};
    }

    FieldBoundaryConditionDivergenceFreeTransverseDirichlet(std::array<value_type, N> const& values,
                                                            DirichletExtrapolation extrapolation)
    {
        for (std::size_t i = 0; i < N; ++i)
            scalar_dirichlet_conditions_[i] = _scalar_dirichlet_bc_type{values[i], extrapolation};
    }

    FieldBoundaryConditionDivergenceFreeTransverseDirichlet(
        FieldBoundaryConditionDivergenceFreeTransverseDirichlet const&) = default;
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet&
    operator=(FieldBoundaryConditionDivergenceFreeTransverseDirichlet const&) = default;
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet(
        FieldBoundaryConditionDivergenceFreeTransverseDirichlet&&) = default;
    FieldBoundaryConditionDivergenceFreeTransverseDirichlet&
    operator=(FieldBoundaryConditionDivergenceFreeTransverseDirichlet&&) = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet;
    }

    void apply(VecFieldT& vecField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               double const time) override
    {
        Direction const direction = getDirection(boundaryLocation);
        Side const side           = getSide(boundaryLocation);
        std::size_t const iNormal = static_cast<std::size_t>(direction);

        auto fields = vecField.components();

        assert(gridLayout.centering(vecField) == gridLayout.centering(tensor_quantity_type::B));

        // handle transverse components with Dirichlet
        for_N<N>([&](auto iTransverse) {
            if (static_cast<std::size_t>(iTransverse) != iNormal)
            {
                field_type& tField = std::get<iTransverse>(fields);
                scalar_dirichlet_conditions_[iTransverse].apply(tField, boundaryLocation,
                                                                localGhostBox, gridLayout, time);
            }
        });

        // set the normal component so the discrete divergence of B is zero, given the transverse
        // ghosts filled above (shared with the transverse-Neumann condition).
        applyDivergenceFreeNormalComponent<dimension>(fields, iNormal, side, gridLayout,
                                                      localGhostBox);
    }

private:
    using _scalar_dirichlet_bc_type = FieldBoundaryConditionDirichlet<field_type, GridLayoutT>;

    std::array<_scalar_dirichlet_bc_type, N> scalar_dirichlet_conditions_;

}; // class FieldBoundaryConditionDivergenceFreeTransverseDirichlet

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_DIVERGENCE_FREE_TRANSVERSE_DIRICHLET_HPP
