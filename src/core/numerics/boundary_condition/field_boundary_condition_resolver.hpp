#ifndef PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_RESOLVER
#define PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_RESOLVER

#include "core/data/tensorfield/tensorfield_traits.hpp"
#include "core/data/vecfield/vecfield_traits.hpp"
#include "core/numerics/boundary_condition/field_antisymmetric_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_dirichlet_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_divergence_free_transverse_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_neumann_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_none_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_symmetric_boundary_condition.hpp"

namespace PHARE::core
{

template<FieldBoundaryConditionType type>
struct FieldBoundaryConditionSelector;

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::None>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldNoneBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::Dirichlet>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldDirichletBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::AntiSymmetric>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldAntiSymmetricBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::Symmetric>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldSymmetricBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::Neumann>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldNeumannBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::DivergenceFreeTransverseNeumann>
{
    // only makes sense for a vector field
    template<IsVecField VecFieldT, typename GridLayoutT>
    using type = FieldDivergenceFreeTransverseNeumannBoundaryCondition<VecFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet>
{
    // only makes sense for a vector field
    template<IsVecField VecFieldT, typename GridLayoutT>
    using type = FieldDivergenceFreeTransverseDirichletBoundaryCondition<VecFieldT, GridLayoutT>;
};

template<FieldBoundaryConditionType type, IsScalarOrTensorField ScalarOrTensorFieldT,
         typename GridLayoutT>
using FieldBoundaryCondition
    = FieldBoundaryConditionSelector<type>::template type<ScalarOrTensorFieldT, GridLayoutT>;

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_RESOLVER
