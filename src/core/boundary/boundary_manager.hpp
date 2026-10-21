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
#include <array>
#include <concepts>
#include <cstdint>
#include <format>
#include <memory>
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
    static std::array<std::string, 3> const dirs{"x", "y", "z"};
    static std::array<std::string, 2> const sides{"lower", "upper"};

    if (!grid.contains("periodicities"))
        return; // no directions declared (minimal
                // test dict): nothing to validate

    for (std::size_t d = 0; d < dimension; ++d)
    {
        if (grid["periodicities"][dirs[d]])
            continue;

        for (auto const& side : sides)
        {
            std::string const loc = std::string{dirs[d]} + side;
            bool const declared
                = grid.contains("domain_boundaries") && grid["domain_boundaries"].contains(loc)
                  && grid["domain_boundaries"][loc].contains("type")
                  && cppdict::get_value(grid["domain_boundaries"][loc], "type", BoundaryType::None)
                         != BoundaryType::None;
            if (!declared)
                throw std::runtime_error(std::format(
                    "BoundaryManager: direction '{0}' is physical but boundary '{1}' has no "
                    "condition declared "
                    "(expected grid/domain_boundaries/{1}/type to be present and not 'none').",
                    dirs[d], loc));
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
    using physical_quantity_type = decltype(GridLayoutT::options.field_options)::Quantity;
    using field_type             = FieldT;
    using boundary_type          = Boundary<GridLayoutT, FieldT>;
    using boundary_factory_type  = BoundaryFactory<GridLayoutT, FieldT>;
    using scalar_quantity_type   = field_type::physical_quantity_type;
    static_assert(std::same_as<scalar_quantity_type, typename physical_quantity_type::Scalar>);
    using vector_field_type     = VecField<field_type, physical_quantity_type>;
    using scalar_condition_type = IFieldBoundaryCondition<field_type, GridLayoutT>;
    using vector_condition_type = IFieldBoundaryCondition<vector_field_type, GridLayoutT>;

    /** @brief Describes how the master boundary is chosen at corner and edges */
    enum class PriorityPolicy : std::uint8_t {
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

        for (auto const& [key, subdict] : dict)
        {
            auto const location   = getBoundaryLocationFromString(key);
            boundaries_[location] = boundary_factory_type::create(
                location, *subdict, scalarQuantities, vectorQuantities, gamma_);
        }
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
    using boundary_map_type = std::unordered_map<BoundaryLocation, std::shared_ptr<boundary_type>>;

    boundary_map_type boundaries_;   //!< List of boundaries mapped by their location.
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
