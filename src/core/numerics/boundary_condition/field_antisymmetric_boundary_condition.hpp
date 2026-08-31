#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_ANTISYMMETRIC_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_ANTISYMMETRIC_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"

namespace PHARE::core
{
/**
 * @brief Anti-symmetric boundary condition for scalar and vector fields.
 *
 * For scalars, imposes a zero value on the boundary (Dirichlet zero).
 * For vectors, imposes zero value on tangential components, Neumann on the normal component.
 *
 * @tparam ScalarOrTensorFieldT Type of the field or tensor field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldAntiSymmetricBoundaryCondition
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static constexpr bool is_scalar        = Super::is_scalar;

    FieldAntiSymmetricBoundaryCondition() = default;

    FieldAntiSymmetricBoundaryCondition(FieldAntiSymmetricBoundaryCondition const&) = default;
    FieldAntiSymmetricBoundaryCondition& operator=(FieldAntiSymmetricBoundaryCondition const&)
        = default;
    FieldAntiSymmetricBoundaryCondition(FieldAntiSymmetricBoundaryCondition&&)            = default;
    FieldAntiSymmetricBoundaryCondition& operator=(FieldAntiSymmetricBoundaryCondition&&) = default;

    virtual ~FieldAntiSymmetricBoundaryCondition() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::AntiSymmetric;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField, BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               double const time) override
    {
        Direction const direction = getDirection(boundaryLocation);

        auto fields = Super::asComponentTuple(scalarOrTensorField);

        for_N<N>([&](auto i) {
            field_type& field = std::get<i>(fields);
            if constexpr (is_scalar)
            {
                scalar_dirichlet_condition_.apply(field, boundaryLocation, localGhostBox,
                                                  gridLayout, time);
            }
            else
            {
                if (static_cast<std::size_t>(i) != static_cast<std::size_t>(direction))
                    scalar_dirichlet_condition_.apply(field, boundaryLocation, localGhostBox,
                                                      gridLayout, time);
                else
                    scalar_neumann_condition_.apply(field, boundaryLocation, localGhostBox,
                                                    gridLayout, time);
            }
        });
    }

private:
    using _scalar_neumann_condition_type = FieldNeumannBoundaryCondition<field_type, GridLayoutT>;
    using _scalar_dirichlet_condition_type
        = FieldDirichletBoundaryCondition<field_type, GridLayoutT>;

    _scalar_neumann_condition_type scalar_neumann_condition_{};
    _scalar_dirichlet_condition_type scalar_dirichlet_condition_{};

}; // class FieldAntiSymmetricBoundaryCondition

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_ANTISYMMETRIC_BOUNDARY_CONDITION_HPP
