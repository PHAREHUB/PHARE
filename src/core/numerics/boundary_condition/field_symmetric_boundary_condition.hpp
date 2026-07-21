#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_SYMMETRIC_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_SYMMETRIC_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"

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
class FieldSymmetricBoundaryCondition
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using tensor_quantity_type = Super::tensor_quantity_type;
    using field_type           = Super::field_type;

    static constexpr size_t dimension = Super::dimension;
    static constexpr size_t N         = Super::N;
    static constexpr bool is_scalar   = Super::is_scalar;

    FieldSymmetricBoundaryCondition() = default;

    FieldSymmetricBoundaryCondition(FieldSymmetricBoundaryCondition const&)            = default;
    FieldSymmetricBoundaryCondition& operator=(FieldSymmetricBoundaryCondition const&) = default;
    FieldSymmetricBoundaryCondition(FieldSymmetricBoundaryCondition&&)                 = default;
    FieldSymmetricBoundaryCondition& operator=(FieldSymmetricBoundaryCondition&&)      = default;

    virtual ~FieldSymmetricBoundaryCondition() = default;

    FieldBoundaryConditionType getType() const override
    {
        return FieldBoundaryConditionType::Symmetric;
    }

    void apply(ScalarOrTensorFieldT& scalarOrTensorField,
               BoundaryLocation const boundaryLocation,
               Box<std::uint32_t, dimension> const& localGhostBox, GridLayoutT const& gridLayout,
               Super::boundary_condition_context_type const& ctx) override
    {
        Direction const direction = getDirection(boundaryLocation);

        auto fields = [&]() {
            if constexpr (is_scalar)
                return std::make_tuple(scalarOrTensorField);
            else
                return scalarOrTensorField.components();
        }();

        for_N<N>([&](auto i) {
            field_type& field = std::get<i>(fields);
            if constexpr (is_scalar)
            {
                scalar_neumann_condition_.apply(field, boundaryLocation, localGhostBox, gridLayout,
                                                ctx);
            }
            else
            {
                if (static_cast<size_t>(i) != static_cast<size_t>(direction))
                    scalar_neumann_condition_.apply(field, boundaryLocation, localGhostBox,
                                                    gridLayout, ctx);
                else
                    scalar_dirichlet_condition_.apply(field, boundaryLocation, localGhostBox,
                                                      gridLayout, ctx);
            }
        });
    }

private:
    using _scalar_neumann_condition_type = FieldNeumannBoundaryCondition<field_type, GridLayoutT>;
    using _scalar_dirichlet_condition_type
        = FieldDirichletBoundaryCondition<field_type, GridLayoutT>;

    _scalar_neumann_condition_type scalar_neumann_condition_{};
    _scalar_dirichlet_condition_type scalar_dirichlet_condition_{};

}; // class FieldSymmetricBoundaryCondition

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_SYMMETRIC_BOUNDARY_CONDITION_HPP
