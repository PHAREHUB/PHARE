#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_MANAGER
#define PHARE_CORE_BOUNDARY_BOUNDARY_MANAGER

#include "core/boundary/boundary.hpp"
#include "core/boundary/boundary_defs.hpp"
#include "core/boundary/boundary_factory.hpp"
#include "core/data/field/field_traits.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/numerics/boundary_condition/field_boundary_condition.hpp"

#include "initializer/data_provider.hpp"

#include <algorithm>
#include <concepts>
#include <memory>
#include <array>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace PHARE::core
{
/**
 * @brief Fail early if the grid declares a physical direction without a boundary
 * condition on each of its two faces.
 *
 * @tparam dimension Number of spatial dimensions.
 * @param grid  The grid sub-dict
 */
template<std::size_t dimension>
inline void validatePhysicalBoundariesDeclared(initializer::PHAREDict const& grid)
{
    static constexpr std::array<char const*, 3> dirs{"x", "y", "z"};
    if (!grid.contains("periodicities"))
        return; // no directions declared (minimal test dict): nothing to validate

    for (std::size_t d = 0; d < dimension; ++d)
    {
        if (grid["periodicities"][dirs[d]].template to<bool>())
            continue;

        for (auto const* side : {"lower", "upper"})
        {
            std::string const loc = std::string{dirs[d]} + side;
            bool const declared
                = grid.contains("boundaries") && grid["boundaries"].contains(loc)
                  && grid["boundaries"][loc].contains("type")
                  && cppdict::get_value(grid["boundaries"][loc], "type", BoundaryType::None)
                         != BoundaryType::None;
            if (!declared)
                throw std::runtime_error("BoundaryManager: direction '" + std::string{dirs[d]}
                                         + "' is physical but boundary '" + loc
                                         + "' has no condition declared (expected grid/boundaries/"
                                         + loc + "/type to be present and not 'none').");
        }
    }
}

/**
 * @brief Manage the lifecycle and retrieval of physical boundary conditions.
 *
 * Store and provide access to boundary condition objects for both
 * scalar and vector fields based on the boundary location and physical quantity.
 *
 * @tparam GridLayoutT The grid layout type.
 */
template<typename GridLayoutT, IsField FieldT>
class BoundaryManager
{
public:
    using physical_quantity_type = GridLayoutT::Quantity;
    using field_type             = FieldT;
    using boundary_type          = Boundary<GridLayoutT, FieldT>;
    using boundary_factory_type  = BoundaryFactory<GridLayoutT, FieldT>;
    using scalar_quantity_type   = field_type::physical_quantity_type;
    static_assert(std::same_as<scalar_quantity_type, typename physical_quantity_type::Scalar>);
    using vector_field_type     = VecField<field_type, physical_quantity_type>;
    using scalar_condition_type = IFieldBoundaryCondition<field_type, GridLayoutT>;
    using vector_condition_type = IFieldBoundaryCondition<vector_field_type, GridLayoutT>;

    /** @brief Describes how the master boundary is chosen at corner and edges */
    enum class PriorityPolicy {
        ByDirection,
        ByBoundaryType,
    };

    BoundaryManager() = delete;

    /**
     * @brief Constructor. Register boundary conditions based on inputfile data.
     * @param dict Configuration dictionary.
     * @param scalar_quantities List of scalar quantities to manage.
     * @param vector_quantities List of vector quantities to manage.
     * @param gamma Heat capacity ratio, required to compute the inflow total energy
     */
    BoundaryManager(PHARE::initializer::PHAREDict const& dict,
                    std::vector<typename physical_quantity_type::Scalar> const& scalarQuantities,
                    std::vector<typename physical_quantity_type::Vector> const& vectorQuantities,
                    double const gamma            = 0.0,
                    PriorityPolicy priorityPolicy = PriorityPolicy::ByDirection)

        : gamma_{gamma}
        , priority_policy_{priorityPolicy}
    {
        if (!dict.isNode())
            return;

        dict.visit(cppdict::visit_all_nodes, [&](std::string const& locationName,
                                                 initializer::PHAREDict::data_t _) {
            /// @todo I don't do anything with the second argument because it cannot be
            /// transformed back into a dict. Maybe add the corresponding constructor to
            /// cppdict, or add the possibility to have a lambda with the second arg
            /// being a dict ?
            BoundaryLocation location = getBoundaryLocationFromString(locationName);
            boundaries_[location]     = boundary_factory_type::create(
                location, dict[locationName], scalarQuantities, vectorQuantities, gamma_);
        });
    }


    /**
     * @brief Retrieve the boundary for a specific location.
     *
     * @param location The location of the desired boundary.
     * @return Non-owning pointer to the matching boundary (owned by this manager), or nullptr if
     *         not found. Returning a raw observer avoids an atomic refcount pair on every lookup
     *         inside the per-boundary-box fill loop.
     */
    boundary_type* getBoundary(BoundaryLocation location) const
    {
        auto it = boundaries_.find(location);
        return (it != boundaries_.end()) ? it->second.get() : nullptr;
    }



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

    _boundary_map_type boundaries_;  //!< List of boundaries mapped by their location.
    double gamma_;                   //!< heat capacity ratio
    PriorityPolicy priority_policy_; //!< How the master boundary is chosen at corners and edges.

    /**
     * @brief Worker function to get the master of an array of 1-codimensional boundary locations,
     * according to the priority policy of the boundary manager.
     *
     * @tparam N Number of elements in the array
     * @param locations Array of boundary locations.
     * @return The location of the master boundary.
     */
    template<std::size_t N>
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
