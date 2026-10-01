#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_SYMMETRIC_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_SYMMETRIC_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_dirichlet.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_neumann.hpp"

namespace PHARE::core
{
/**
 * @brief Symmetric boundary condition for scalar and vector fields.
 *
 * For scalars, imposes a null derivative along the normal (Neumann).
 * For vectors, imposes Neumann on tangential components, zero value on the normal component.
 *
 * @tparam ScalarOrTensorFieldT Type of the field or tensor field.
 * @tparam GridLayoutT Grid layout configuration.
 *
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldBoundaryConditionSymmetric
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;

    static constexpr std::size_t dimension = Super::dimension;
    static constexpr std::size_t N         = Super::N;
    static constexpr bool is_scalar        = Super::is_scalar;

    FieldBoundaryConditionSymmetric() = default;

    FieldBoundaryConditionSymmetric(FieldBoundaryConditionSymmetric const&)            = default;
    FieldBoundaryConditionSymmetric& operator=(FieldBoundaryConditionSymmetric const&) = default;
    FieldBoundaryConditionSymmetric(FieldBoundaryConditionSymmetric&&)                 = default;
    FieldBoundaryConditionSymmetric& operator=(FieldBoundaryConditionSymmetric&&)      = default;

    virtual ~FieldBoundaryConditionSymmetric() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Symmetric;
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
                scalar_neumann_condition_.apply(field, boundaryLocation, localGhostBox, gridLayout,
                                                time);
            }
            else
            {
                if (static_cast<std::size_t>(i) != static_cast<std::size_t>(direction))
                    scalar_neumann_condition_.apply(field, boundaryLocation, localGhostBox,
                                                    gridLayout, time);
                else
                    scalar_dirichlet_condition_.apply(field, boundaryLocation, localGhostBox,
                                                      gridLayout, time);
            }
        });
    }

private:
    using _scalar_neumann_condition_type = FieldBoundaryConditionNeumann<field_type, GridLayoutT>;
    using _scalar_dirichlet_condition_type
        = FieldBoundaryConditionDirichlet<field_type, GridLayoutT>;

    _scalar_neumann_condition_type scalar_neumann_condition_{};
    _scalar_dirichlet_condition_type scalar_dirichlet_condition_{};

}; // class FieldBoundaryConditionSymmetric

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_SYMMETRIC_HPP
