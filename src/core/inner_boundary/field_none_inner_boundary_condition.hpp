#ifndef PHARE_CORE_INNER_BOUNDARY_FIELD_NONE_INNER_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_INNER_BOUNDARY_FIELD_NONE_INNER_BOUNDARY_CONDITION_HPP

#include "core/inner_boundary/field_inner_boundary_condition.hpp"

namespace PHARE::core
{
template<typename ScalarOrTensorFieldT, typename GridLayoutT, typename PhysicalStateT>
class FieldNoneInnerBoundaryCondition
    : public FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>
{
public:
    using Super = FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    using inner_boundary_mesh_data_type = Super::inner_boundary_mesh_data_type;
    using context_type                  = Super::context_type;

    FieldNoneInnerBoundaryCondition() = default;

    FieldInnerBoundaryConditionType getType() const override
    {
        return FieldInnerBoundaryConditionType::None;
    }

    void apply(ScalarOrTensorFieldT& /*scalarOrTensorField*/, GridLayoutT const& /*layout*/,
               inner_boundary_mesh_data_type const& /*boundaryMeshData*/,
               context_type const& /*ctx*/) override
    {
    }
};

} // namespace PHARE::core
#endif // PHARE_CORE_INNER_BOUNDARY_FIELD_NONE_INNER_BOUNDARY_CONDITION_HPP
