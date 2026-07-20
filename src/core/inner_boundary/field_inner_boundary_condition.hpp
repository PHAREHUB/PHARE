#ifndef PHARE_CORE_INNER_BOUNDARY_FIELD_INNER_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_INNER_BOUNDARY_FIELD_INNER_BOUNDARY_CONDITION_HPP

#include "core/data/field/field_traits.hpp"
#include "core/data/tensorfield/tensorfield_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/inner_boundary/inner_bc_context.hpp"
#include "core/inner_boundary/inner_boundary_geometry.hpp"
#include "core/inner_boundary/inner_boundary_mesh_data.hpp"
#include "core/numerics/interpolator/field_at_point.hpp"

namespace PHARE::core
{

/**
 * @brief Taxonomy of supported inner-boundary condition kinds for field quantities.
 *
 * - **Dirichlet** — the field value is prescribed on the boundary (e.g. impose B·n = 0).
 * - **Neumann** — the normal derivative of the field is prescribed on the boundary.
 * - **Symmetric** — the field is mirrored across the boundary (even reflection).
 * - **Antisymmetric** — the field is sign-flipped across the boundary (odd reflection).
 * - **AdaptiveDirichletOrNeumann** — Dirichlet or Neumann selected per ghost element from the
 *   sign of a target vector quantity projected on the boundary normal (characteristic switch).
 */
enum class FieldInnerBoundaryConditionType {
    None,
    Dirichlet,
    Neumann,
    Symmetric,
    Antisymmetric,
    AdaptiveDirichletOrNeumann,
    IonosphericConvectionMomentum
};

/**
 * @brief Abstract base class for applying a boundary condition to a scalar or tensor field
 *        at an inner boundary.
 *
 * Concrete subclasses implement `operator()` to enforce a specific condition
 * (Dirichlet, Neumann, Symmetric, or Antisymmetric) by writing corrected values into
 * the ghost mesh elements identified by the associated @ref InnerBoundaryMeshData.
 *
 *
 * @tparam ScalarOrTensorFieldT Either a type satisfying `IsField` (scalar) or `IsTensorField`.
 * @tparam GridLayoutT          Grid layout providing geometry and centering information.
 */
template<IsScalarOrTensorField ScalarOrTensorFieldT, typename GridLayoutT, typename PhysicalStateT>
class FieldInnerBoundaryCondition
{
public:
    static constexpr bool is_scalar   = IsField<ScalarOrTensorFieldT>;
    static constexpr size_t dimension = GridLayoutT::dimension;
    static constexpr size_t N = NumberOfComponentsSelector<ScalarOrTensorFieldT, is_scalar>::value;

    using This = FieldInnerBoundaryCondition<ScalarOrTensorFieldT, GridLayoutT, PhysicalStateT>;

    using inner_boundary_type    = InnerBoundaryGeometry<dimension>;
    using field_type             = FieldTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;
    using physical_quantity_type = GridLayoutT::Quantity;
    using vecfield_type          = VecField<field_type, physical_quantity_type>;
    using tensor_quantity_type
        = PhysicalQuantityTypeSelector<ScalarOrTensorFieldT, is_scalar>::type;
    using inner_boundary_mesh_data_type = InnerBoundaryMeshData<dimension, physical_quantity_type>;
    using ghost_elem_data_type          = inner_boundary_mesh_data_type::ghost_elem_data_type;
    using state_type                    = PhysicalStateT;
    using context_type                  = InnerBCContext<PhysicalStateT>;
    using interpolator_type
        = FieldAtPoint<dimension,
                       /*interpOrder=*/1>; ///< Point interpolator with hard-coded linear
                                           ///< order to guaranty that interpolation
                                           ///< support does not cross the boundary. Going-higher
                                           ///< order would require directional interpolation.

    FieldInnerBoundaryCondition()          = default;
    virtual ~FieldInnerBoundaryCondition() = default;

    FieldInnerBoundaryCondition(This const&) = delete;
    FieldInnerBoundaryCondition(This&&)      = delete;
    This& operator=(This const&)             = delete;
    This& operator=(This&&)                  = delete;

    /**
     * @brief Returns the kind of boundary condition implemented by this object.
     */
    virtual FieldInnerBoundaryConditionType getType() const = 0;

    /**
     * @brief Application priority used by InnerBoundaryManager::applyToMoments to order BCs.
     *
     * Lower runs first. Most BCs are independent and keep the default (0). A BC that
     * depends on other quantities' ghost cells already being filled (e.g. total energy
     * reconstructed from pressure + ghost rho/rhoV) returns a higher value so it runs last.
     */
    virtual int priority() const { return 0; }

    /**
     * @brief Whether this BC writes a value into ghost cells that have no interpolable
     * fluid-side sample ("non-settable" ghosts).
     *
     * Most BCs need to interpolate at the symmetric point and so leave such ghosts untouched
     * (default false). A 0th-order (constant) Dirichlet needs no interpolation and fills every
     * ghost; BCs that reconstruct from another quantity's ghosts (e.g. total energy from
     * pressure) query this to decide whether those ghosts are safe to reconstruct.
     */
    virtual bool fillsNonInterpolableGhosts() const { return false; }

    /**
     * @brief Apply the boundary condition to @p scalarOrTensorField.
     *
     * Implementations fill ghost-cell values in @p scalarOrTensorField so that the
     * chosen BC is satisfied at @p boundary.
     *
     * @param scalarOrTensorField Field whose ghost cells are to be updated in place.
     * @param boundary            Inner boundary providing geometry (signed-distance, normal).
     * @param layout              Grid layout used to convert between local and physical indices.
     * @param time                Current simulation time (needed for time-dependent Dirichlet BCs).
     */
    virtual void apply(ScalarOrTensorFieldT& scalarOrTensorField, GridLayoutT const& layout,
                       inner_boundary_mesh_data_type const& boundaryMeshData,
                       context_type const& ctx)
        = 0;

protected:
    interpolator_type interpolator_;
};

} // namespace PHARE::core


#endif // PHARE_CORE_INNER_BOUNDARY_FIELD_INNER_BOUNDARY_CONDITION_HPP
