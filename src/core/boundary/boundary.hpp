#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_HPP
#define PHARE_CORE_BOUNDARY_BOUNDARY_HPP

#include "core/boundary/boundary_defs.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition_factory.hpp"

#include <concepts>
#include <memory>
#include <unordered_map>
#include <utility>

namespace PHARE::core
{
/**
 * @brief A Boundary is associated to one of the physical boundary (XLower, XUpper) and manages the
 * collection of boundary conditions associated with each physical quantities that requires one.
 * This class is not polymorphic itself, but the different type of boundaries are obtained thanks to
 * the polymorphism of the boundary conditions. Which runtime type to choose for a boundary
 * condition applied to a physical quantity is controled by the @cBoundaryFactory following the
 * desired @c BoundaryType.
 *
 * @tparam PhysicalQuantityT The model category of physical quantities (MHDQuantity or
 * HybridQuantity).
 * @tparam FieldT The type for scalar fields.
 * @tparam GridLayoutT The type for the grid layout.
 */
template<typename PhysicalQuantityT, IsField FieldT, IsGridLayout GridLayoutT>
class Boundary
{
public:
    using This                 = Boundary<PhysicalQuantityT, FieldT, GridLayoutT>;
    using scalar_quantity_type = FieldT::physical_quantity_type;
    static_assert(std::same_as<scalar_quantity_type, typename PhysicalQuantityT::Scalar>);
    using vector_quantity_type        = PhysicalQuantityT::Vector;
    using vector_field_type           = VecField<FieldT, PhysicalQuantityT>;
    using scalar_field_condition_type = IFieldBoundaryCondition<FieldT, GridLayoutT>;
    using vector_field_condition_type = IFieldBoundaryCondition<vector_field_type, GridLayoutT>;

    Boundary() = delete;
    Boundary(BoundaryType type, BoundaryLocation location)
        : type_{type}
        , location_{location} {};
    ~Boundary() = default;

    inline BoundaryType getType() const { return type_; };
    inline BoundaryLocation getLocation() const { return location_; };

    /**
     * @brief Retrieve the registered field boundary condition corresponding to a physical quantity.
     *
     * @tparam TensorPhysicalQuantityT Type of the physical quantity, expected to be either @c
     * PhysicalQuantity::Scalar, or @c PhysicalQuantityT::Vector.
     *
     * @param quantity The physical quantity whose field boundary condition is wanted
     * @return A shared pointer to the boundary condition if a one has been previously registered
     *         for the physical quantity, nullptr otherwise.
     */
    template<typename TensorPhysicalQuantityT>
    auto getFieldCondition(TensorPhysicalQuantityT quantity) const
    {
        if constexpr (std::same_as<TensorPhysicalQuantityT, typename PhysicalQuantityT::Scalar>)
        {
            auto it = scalar_field_conditions_.find(quantity);
            return (it != scalar_field_conditions_.end()) ? it->second : nullptr;
        }
        else if constexpr (std::same_as<TensorPhysicalQuantityT,
                                        typename PhysicalQuantityT::Vector>)
        {
            auto it = vector_field_conditions_.find(quantity);
            return (it != vector_field_conditions_.end()) ? it->second : nullptr;
        }
        else
        {
            static_assert(dependant_false_<TensorPhysicalQuantityT>,
                          "Tensoriality of the physical quantity not supported.");
        }
    }

    /**
     * @brief Register a field boundary condition for a quantity to the boundary.
     *
     * @tparam type The corresponding value of the @c BoundaryType enum corresponding to the desired
     * boundary type.
     * @tparam TensorPhysicalQuantityT Type of the physical quantity, expected to be either @c
     * PhysicalQuantity::Scalar, or @c PhysicalQuantityT::Vector.
     * @tparam Args Types of the arguments for the FieldConditionT constructor.
     *
     * @param quantity The physical quantity (scalar or vector) to which the field condition should
     * apply.
     * @param args The arguments for the field BC constructor, passed by perfect forwarding to the
     * field BC factory.
     */
    template<FieldBoundaryConditionType type, typename TensorPhysicalQuantityT, typename... Args>
    void registerFieldCondition(TensorPhysicalQuantityT quantity, Args&&... args)
    {
        if constexpr (std::same_as<TensorPhysicalQuantityT, scalar_quantity_type>)
        {
            scalar_field_conditions_[quantity]
                = FieldBoundaryConditionFactory::create<type, FieldT, GridLayoutT>(
                    std::forward<Args>(args)...);
        }
        else if constexpr (std::same_as<TensorPhysicalQuantityT, vector_quantity_type>)
        {
            vector_field_conditions_[quantity]
                = FieldBoundaryConditionFactory::create<type, vector_field_type, GridLayoutT>(
                    std::forward<Args>(args)...);
        }
        else
        {
            static_assert(dependant_false_<TensorPhysicalQuantityT>,
                          "Tensoriality of the physical quantity not supported.");
        }
    }

    /**
     * @brief Define comparison of boundaries based on the enum @c BoundaryType .
     */
    std::strong_ordering operator<=>(This const& other) const
    {
        return this->getType() <=> other.getType();
    }

private:
    using _scalar_field_condition_map_type
        = std::unordered_map<scalar_quantity_type, std::shared_ptr<scalar_field_condition_type>>;
    using _vector_field_condition_map_type
        = std::unordered_map<vector_quantity_type, std::shared_ptr<vector_field_condition_type>>;

    /** Utility to make compilation fail in certain conditions. */
    template<typename...>
    static constexpr bool dependant_false_ = false;

    /** The type of the boundary (open, inflow, reflective ...) */
    BoundaryType type_;
    /** The location of the boundary (XLower, XUpper, ...) */
    BoundaryLocation location_;
    /** The list of registered scalar field boundary conditions on the boundary */
    _scalar_field_condition_map_type scalar_field_conditions_;
    /** The list of registered vector field boundary conditions on the boundary */
    _vector_field_condition_map_type vector_field_conditions_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_BOUNDARY_HPP
