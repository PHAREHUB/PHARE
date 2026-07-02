#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NONE_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NONE_BOUNDARY_CONDITION_HPP

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
class FieldNoneBoundaryCondition
    : public IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>
{
public:
    using Super                    = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    static constexpr size_t dimension = Super::dimension;

    FieldNoneBoundaryCondition() = default;

    FieldNoneBoundaryCondition(FieldNoneBoundaryCondition const&)            = default;
    FieldNoneBoundaryCondition& operator=(FieldNoneBoundaryCondition const&) = default;
    FieldNoneBoundaryCondition(FieldNoneBoundaryCondition&&)                 = default;
    FieldNoneBoundaryCondition& operator=(FieldNoneBoundaryCondition&&)      = default;

    virtual ~FieldNoneBoundaryCondition() = default;

    FieldBoundaryConditionType getType() const override { return FieldBoundaryConditionType::None; }

    void apply(ScalarOrTensorFieldT& /*scalarOrTensorField*/,
               BoundaryLocation const /*boundaryLocation*/,
               Box<std::uint32_t, dimension> const& /*localGhostBox*/,
               GridLayoutT const& /*gridLayout*/,
               [[maybe_unused]] Super::boundary_condition_context_type const& /*ctx*/) override
    {
    }
}; // class FieldNoneBoundaryCondition

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_NONE_BOUNDARY_CONDITION_HPP
