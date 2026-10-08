#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP
#define PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP

#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace PHARE::core
{
/** @brief Physical behavior of a boundary. */
enum class BoundaryType { None, Reflective, SuperMagnetofastInflow, Open };

/**
 * @brief Possible locations of 1-codimensional boundary (a face in 3D, an edge in 2D, an extremity
 * in 1D).
 */
enum class BoundaryLocation {
    XLower = 0,
    XUpper = 1,
    YLower = 2,
    YUpper = 3,
    ZLower = 4,
    ZUpper = 5
};

/**
 * @brief Return the side of a boundary location.
 * @param boundaryLoc The boundary location.
 * @return The boundary side.
 */
constexpr Side getSide(BoundaryLocation boundaryLoc)
{
    switch (boundaryLoc)
    {
        case BoundaryLocation::XLower:
        case BoundaryLocation::YLower:
        case BoundaryLocation::ZLower: return Side::Lower; break;

        case BoundaryLocation::XUpper:
        case BoundaryLocation::YUpper:
        case BoundaryLocation::ZUpper: return Side::Upper; break;

        default: throw std::runtime_error("Invalid BoundaryLocation.");
    }
};

/** @brief Return the direction of a boundary location.
 * @param boundaryLoc The boundary location.
 * @return The boundary direction.
 */
constexpr Direction getDirection(BoundaryLocation boundaryLoc)
{
    switch (boundaryLoc)
    {
        case BoundaryLocation::XLower:
        case BoundaryLocation::XUpper: return Direction::X; break;

        case BoundaryLocation::YLower:
        case BoundaryLocation::YUpper: return Direction::Y; break;

        case BoundaryLocation::ZLower:
        case BoundaryLocation::ZUpper: return Direction::Z; break;

        default: throw std::runtime_error("Invalid BoundaryLocation.");
    }
};

/** @brief Possible locations of a 2-codimensional boundary (an edge in 3D, a corner in 2D) */
enum class Codim2BoundaryLocation : std::uint16_t {
    XLower_YLower = 0,
    XUpper_YLower = 1,
    XLower_YUpper = 2,
    XUpper_YUpper = 3,
    XLower_ZLower = 4,
    XUpper_ZLower = 5,
    XLower_ZUpper = 6,
    XUpper_ZUpper = 7,
    YLower_ZLower = 8,
    YUpper_ZLower = 9,
    YLower_ZUpper = 10,
    YUpper_ZUpper = 11
};

/**
 * @brief Return the location of the two (1-codimensional) boundaries adjacent to a 2-codimensional
 * boundary.
 * @param The location of the 2-codimensional boundary.
 * @return An array containing the two locations of the adjacent boundaries.
 */
constexpr std::array<BoundaryLocation, 2>
getAdjacentBoundaryLocations(Codim2BoundaryLocation const location)
{
    using enum BoundaryLocation;
    std::array<std::array<BoundaryLocation, 2>, 12> constexpr adjacents{{
        {XLower, YLower},
        {XUpper, YLower},
        {XLower, YUpper},
        {XUpper, YUpper},
        {XLower, ZLower},
        {XUpper, ZLower},
        {XLower, ZUpper},
        {XUpper, ZUpper},
        {YLower, ZLower},
        {YUpper, ZLower},
        {YLower, ZUpper},
        {YUpper, ZUpper},
    }};

    auto const idx = static_cast<std::underlying_type_t<Codim2BoundaryLocation>>(location);

    if (idx >= adjacents.size())
        throw std::runtime_error("Invalid adjacent boundary location index.");

    return adjacents[idx];
}

/** @brief Possible locations of a 3-codimensional boundary (a corner in 3D) */
enum class Codim3BoundaryLocation : std::uint16_t {
    XLower_YLower_ZLower = 0,
    XUpper_YLower_ZLower = 1,
    XLower_YUpper_ZLower = 2,
    XUpper_YUpper_ZLower = 3,
    XLower_YLower_ZUpper = 4,
    XUpper_YLower_ZUpper = 5,
    XLower_YUpper_ZUpper = 6,
    XUpper_YUpper_ZUpper = 7
};

/**
 * @brief Return the location of the three (1-codimensional) boundaries adjacent to a
 * 3-codimensional boundary.
 * @param The location of the 3-codimensional boundary.
 * @return An array containing the three locations of the adjacent boundaries.
 */
constexpr std::array<BoundaryLocation, 3>
getAdjacentBoundaryLocations(Codim3BoundaryLocation const location)
{
    using enum BoundaryLocation;
    std::array<std::array<BoundaryLocation, 3>, 8> constexpr adjacents{{
        {XLower, YLower, ZLower},
        {XUpper, YLower, ZLower},
        {XLower, YUpper, ZLower},
        {XUpper, YUpper, ZLower},
        {XLower, YLower, ZUpper},
        {XUpper, YLower, ZUpper},
        {XLower, YUpper, ZUpper},
        {XUpper, YUpper, ZUpper},
    }};

    auto const idx = static_cast<std::underlying_type_t<Codim3BoundaryLocation>>(location);

    if (idx >= adjacents.size())
        throw std::runtime_error("Invalid adjacent boundary location index.");

    return adjacents[idx];
}

/**
 * @brief Get the BoundaryLocation from input keyword, and throw an error if the keyword does not
 * correspond to any known boundary location.
 */
inline BoundaryLocation getBoundaryLocationFromString(std::string const& name)
{
    static std::unordered_map<std::string, BoundaryLocation> const typeMap_ = {
        {"xlower", BoundaryLocation::XLower}, {"xupper", BoundaryLocation::XUpper},
        {"ylower", BoundaryLocation::YLower}, {"yupper", BoundaryLocation::YUpper},
        {"zlower", BoundaryLocation::ZLower}, {"zupper", BoundaryLocation::ZUpper},
    };

    auto it = typeMap_.find(name);
    if (it == typeMap_.end())
        throw std::runtime_error("Wrong boundary location name = " + name);
    return it->second;
}

/**
 * @brief Meta utilities to retrieve the enum type of boundary location depending on the
 * codimension.
 * @tparam N Codimension value.
 */
template<std::size_t N>
using CodimNBoundaryLocation = std::tuple_element_t<
    N - 1, std::tuple<BoundaryLocation, Codim2BoundaryLocation, Codim3BoundaryLocation>>;


} // namespace PHARE::core

#endif /* PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP */
