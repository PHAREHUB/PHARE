#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_MANAGER
#define PHARE_CORE_BOUNDARY_BOUNDARY_MANAGER

#include "core/boundary/boundary.hpp"
#include "core/boundary/boundary_defs.hpp"
#include "core/boundary/boundary_factory.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/grid/gridlayout_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include "core/numerics/thermo/thermo.hpp"
#include "initializer/data_provider.hpp"

#include <algorithm>
#include <concepts>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace PHARE::core
{
/**
 * @brief Manage the lifecycle and retrieval of physical boundary conditions.
 *
 * Store and provide access to boundary condition objects for both
 * scalar and vector fields based on the boundary location and physical quantity.
 *
 * @tparam PhysicalQuantityT Type defining scalar and vector quantities (MHDQuantity or
 * HybridQuantity).
 * @tparam FieldT The scalar field type.
 * @tparam GridLayoutT The grid layout type.
 */
template<typename PhysicalQuantityT, IsField FieldT, IsGridLayout GridLayoutT>
class BoundaryManager
{
public:
    using boundary_type          = Boundary<PhysicalQuantityT, FieldT, GridLayoutT>;
    using boundary_factory_type  = BoundaryFactory<PhysicalQuantityT, FieldT, GridLayoutT>;
    using physical_quantity_type = PhysicalQuantityT;
    using scalar_quantity_type   = FieldT::physical_quantity_type;
    static_assert(std::same_as<scalar_quantity_type, typename physical_quantity_type::Scalar>);
    using vector_field_type     = VecField<FieldT, PhysicalQuantityT>;
    using scalar_condition_type = IFieldBoundaryCondition<FieldT, GridLayoutT>;
    using vector_condition_type = IFieldBoundaryCondition<vector_field_type, GridLayoutT>;

    BoundaryManager() = delete;

    /**
     * @brief Constructor. Register boundary conditions based on inputfile data.
     * @param dict Configuration dictionary.
     * @param scalar_quantities List of scalar quantities to manage.
     * @param vector_quantities List of vector quantities to manage.
     * @param thermo Optional thermodynamic model, required for EOS-dependent boundary conditions
     *               (e.g. inflow). Pass nullptr (default) for models that do not use an EOS.
     */
    BoundaryManager(PHARE::initializer::PHAREDict const& dict,
                    std::vector<typename PhysicalQuantityT::Scalar> const& scalarQuantities,
                    std::vector<typename PhysicalQuantityT::Vector> const& vectorQuantities,
                    std::shared_ptr<Thermo> thermo = nullptr)
        : thermo_{std::move(thermo)}
    {
        dict.visit(cppdict::visit_all_nodes,
                   [&](std::string const& locationName, initializer::PHAREDict::data_t _) {
                       /// @todo I don't do anything with the second argument because it cannot be
                       /// transformed back into a dict. Maybe add the corresponding constructor to
                       /// cppdict, or add the possibility to have a lambda with the second arg
                       /// being a dict ?
                       BoundaryLocation location = getBoundaryLocationFromString(locationName);
                       boundaries_[location]     = boundary_factory_type::create(
                           location, dict[locationName], scalarQuantities, vectorQuantities,
                           thermo_);
                   });

        /// @todo If this mode stays in the code it should be read from the input dict.
        priority_policy_ = PriorityPolicy::ByDirection;
    }


    /**
     * @brief Retrieve the boundary for a specific location.
     *
     * @param location The location of the desired boundary.
     * @return Shared pointer to the matching boundary, or nullptr if not found.
     *
     */
    std::shared_ptr<boundary_type> getBoundary(BoundaryLocation location) const
    {
        auto it = boundaries_.find(location);
        return (it != boundaries_.end()) ? it->second : nullptr;
    }

    /** @brief Describes how the master boundary is chosen at corner and edges */
    enum class PriorityPolicy {
        ByDirection,
        ByBoundaryType,
    };

    void setPriorityPolicy(PriorityPolicy policy) { priority_policy_ = policy; }

    /** @brief Gets the master 1-codimensional boundary for any given N-codimensional boundary,
     * following the priority policy of the boundary manager.
     *
     * @note If @p location corresponds itself to a 1-codim boundary, then it returns the same
     * @p location.
     *
     * @tparam CodimNBoundaryLocationT Type of boundary location.
     * @param location The location of the boundary where we want to determine which is the master
     * boundary.
     * @return The location of the master boundary.
     */
    template<typename CodimNBoundaryLocationT>
    BoundaryLocation getMasterBoundaryLocation(CodimNBoundaryLocationT location) const
    {
        if constexpr (std::same_as<CodimNBoundaryLocationT, BoundaryLocation>)
        {
            return location;
        }
        else
        {
            return selectMasterBoundaryInArray_(getAdjacentBoundaryLocations(location));
        }
    }

private:
    using _boundary_map_type = std::unordered_map<BoundaryLocation, std::shared_ptr<boundary_type>>;

    /** @brief Utility struct to group scalar and vector quantities together */
    struct SimulationMenu
    {
        std::vector<typename PhysicalQuantityT::Scalar> const& scalars;
        std::vector<typename PhysicalQuantityT::Vector> const& vectors;
    };


    _boundary_map_type boundaries_;  //!< List of boundaries mapped by their location.
    PriorityPolicy priority_policy_; //!< How the master boundary is chosen at corners and edges.
    std::shared_ptr<Thermo> thermo_; //!< EOS object, nullptr for non-MHD models.

    /**
     * @brief Worker function to get the master of an array of 1-codimensional boundary locations,
     * according to the priority policy of the boundary manager.
     *
     * @tparam N Number of elements in the array
     * @param locations Array of boundary locations.
     * @return The location of the master boundary.
     */
    template<size_t N>
    BoundaryLocation selectMasterBoundaryInArray_(std::array<BoundaryLocation, N> locations) const
    {
        switch (priority_policy_)
        {
            case PriorityPolicy::ByDirection: {
                auto it = std::ranges::max_element(locations, {}, getDirection);
                return *it;
            }

            case PriorityPolicy::ByBoundaryType: {
                auto it = std::ranges::max_element(locations, {}, [&](auto location) {
                    if (auto boundaryPtr = getBoundary(location); boundaryPtr)
                        return boundaryPtr->getType();
                    else
                        throw std::runtime_error("Pointer to boundary is null.");
                });
                return *it;
            }

            default: throw std::runtime_error("Non-existing priority mode for boundaries.");
        }
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_BOUNDARY_BOUNDARY_MANAGER
