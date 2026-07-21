#ifndef PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/patch_field_accessor.hpp"
#include "core/data/tensorfield/tensorfield_traits.hpp"
#include "core/numerics/boundary_condition/boundary_condition_context.hpp"
#include "core/utilities/box/box.hpp"

namespace PHARE::core
{

/**
 * @brief Supported types of field boundary conditions.
 *
 */
enum class FieldBoundaryConditionType : int {
    None,
    Dirichlet,
    AntiSymmetric,
    Symmetric,
    Neumann,
    DivergenceFreeTransverseNeumann,
    DivergenceFreeTransverseDirichlet,
    TotalEnergyFromPressure
};


/**
 * @brief Interface for applying boundary conditions to scalar or tensor fields.
 *
 * Concrete field boundary conditions are provided by implementating this interface.
 *
 * @tparam ScalarOrTensorFieldT The type of the scalarOrTensorField (must satisfy IsField or
 * IsTensorField).
 * @tparam GridLayoutT The grid layout type (must satisfy IsGridLayout).
 *
 */
template<typename ScalarOrTensorFieldT, IsGridLayout GridLayoutT>
    requires(IsField<ScalarOrTensorFieldT> || IsTensorField<ScalarOrTensorFieldT>)
class IFieldBoundaryCondition
{
public:
    static constexpr bool is_scalar   = IsField<ScalarOrTensorFieldT>;
    static constexpr size_t dimension = GridLayoutT::dimension;
    static constexpr size_t N = NumberOfComponentsSelector<ScalarOrTensorFieldT, is_scalar>::value;

    using This                   = IFieldBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT>;
    using physical_quantity_type = typename GridLayoutT::Quantity;
    using tensor_quantity_type
        = PhysicalQuantityTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;
    using field_type = FieldTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;
    using patch_field_accessor_type
        = IPatchFieldAccessor<field_type, typename GridLayoutT::Quantity>;
    using boundary_condition_context_type
        = BoundaryConditionContext<field_type, typename GridLayoutT::Quantity>;


    /** @brief Return the type of the boundary condition. */
    virtual FieldBoundaryConditionType getType() const = 0;

    virtual ~IFieldBoundaryCondition() = default;


    /**
     * @brief Enforce the boundary condition on the provided scalar/tensor @p scalarOrTensorField,
     * by filling accordingly the ghost cells contained in the local box @p localGhostBox, at the
     * physical time carried by @p ctx, and considering that the boundary is located at
     * @p boundaryLocation.
     *
     * @param scalarOrTensorField The scalar or tensor to which we apply the boundary condition.
     * @param boundaryLocation The location of the physical boundary.
     * @param localGhostBox The box containing the ghost cells/nodes to fill.
     * @param gridLayout The grid layout.
     * @param ctx Bundle of context data: accessors to the current and previous substage states,
     *            simulation time, and substage time step. State-aware BCs read the previous-state
     *            accessor `ctx.accessor_old` and write into ghost cells reachable via
     *            `ctx.accessor_new`; simple BCs only need `ctx.accessor_new` and `ctx.time`.
     */
    virtual void apply(ScalarOrTensorFieldT& scalarOrTensorField,
                       BoundaryLocation const boundaryLocation,
                       Box<std::uint32_t, dimension> const& localGhostBox,
                       GridLayoutT const& gridLayout, boundary_condition_context_type const& ctx)
        = 0;
};

} // namespace PHARE::core
#endif // PHARE_CORE_NUMERICS_BOUNDARY_CONDITION_FIELD_BOUNDARY_CONDITION_HPP
