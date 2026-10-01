#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_NONE_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_NONE_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include <cstddef>

namespace PHARE::core
{
/**
 * @brief 'None' boundary condition for scalar and vector fields.
 *
 * @tparam ScalarOrTensorFieldT Type of the field or tensor field.
 * @tparam GridLayoutT Grid layout configuration.
 */
template<typename ScalarOrTensorFieldT, typename GridLayoutT>
class FieldBoundaryConditionNone : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    static constexpr std::size_t dimension = Super::dimension;

    FieldBoundaryConditionNone() = default;

    FieldBoundaryConditionNone(FieldBoundaryConditionNone const&)            = default;
    FieldBoundaryConditionNone& operator=(FieldBoundaryConditionNone const&) = default;
    FieldBoundaryConditionNone(FieldBoundaryConditionNone&&)                 = default;
    FieldBoundaryConditionNone& operator=(FieldBoundaryConditionNone&&)      = default;

    virtual ~FieldBoundaryConditionNone() = default;

    FieldBoundaryConditionType getType() const override { return FieldBoundaryConditionType::None; }

    void apply(ScalarOrTensorFieldT& /*scalarOrTensorField*/,
               BoundaryLocation const /*boundaryLocation*/,
               Box<std::uint32_t, dimension> const& /*localGhostBox*/,
               GridLayoutT const& /*gridLayout*/, [[maybe_unused]] double const /*time*/) override
    {
    }
}; // class FieldBoundaryConditionNone

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_NONE_HPP
