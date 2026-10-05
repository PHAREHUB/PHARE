#ifndef PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_RESOLVER
#define PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_RESOLVER

#include "core/data/tensorfield/tensorfield_traits.hpp"
#include "core/data/vecfield/vecfield_traits.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_antisymmetric.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_dirichlet.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_divergence_free_transverse_dirichlet.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_divergence_free_transverse_neumann.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_neumann.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_none.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_symmetric.hpp"

namespace PHARE::core
{

template<FieldBoundaryConditionType type>
struct FieldBoundaryConditionSelector;

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::None>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionNone<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::Dirichlet>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionDirichlet<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::AntiSymmetric>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionAntiSymmetric<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::Symmetric>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionSymmetric<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::Neumann>
{
    template<typename ScalarOrTensorFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionNeumann<ScalarOrTensorFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::DivergenceFreeTransverseNeumann>
{
    // only makes sense for a vector field
    template<IsVecField VecFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionDivergenceFreeTransverseNeumann<VecFieldT, GridLayoutT>;
};

template<>
struct FieldBoundaryConditionSelector<FieldBoundaryConditionType::DivergenceFreeTransverseDirichlet>
{
    // only makes sense for a vector field
    template<IsVecField VecFieldT, typename GridLayoutT>
    using type = FieldBoundaryConditionDivergenceFreeTransverseDirichlet<VecFieldT, GridLayoutT>;
};

template<FieldBoundaryConditionType type, IsScalarOrTensorField ScalarOrTensorFieldT,
         typename GridLayoutT>
using FieldBoundaryCondition
    = FieldBoundaryConditionSelector<type>::template type<ScalarOrTensorFieldT, GridLayoutT>;

} // namespace PHARE::core

#endif // PHARE_CORE_NUMERICS_FIELD_BOUNDARY_CONDITION_RESOLVER
