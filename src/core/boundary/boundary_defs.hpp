#ifndef PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP
#define PHARE_CORE_BOUNDARY_BOUNDARY_DEFS_HPP

#include "core/data/grid/gridlayoutdefs.hpp"

#include <unordered_map>

namespace PHARE::core
{
/** @brief Physical behavior of a boundary. */
enum class BoundaryType {
    None,
    Reflective,
    SuperMagnetofastInflow,
    SuperMagnetofastOutflow,
    Open,
    FreePressureInflow,
    FixedPressureOutflow
};

/** @brief Possible codimension of a boundary. */
enum class BoundaryCodim { One = 1, Two = 2, Three = 3 };

//@{
//! @name Definitions for boundary array sizes in 1d, 2d, or 3d:
int const NUM_1D_NODES = 2;

int const NUM_2D_EDGES = 4;
int const NUM_2D_NODES = 4;

int const NUM_3D_FACES = 6;
int const NUM_3D_EDGES = 12;
int const NUM_3D_NODES = 8;
//@}

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
enum class Codim2BoundaryLocation {
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
getAdjacentBoundaryLocations(Codim2BoundaryLocation location)
{
    switch (location)
    {
        // X-Y Edges
        case Codim2BoundaryLocation::XLower_YLower:
            return {BoundaryLocation::XLower, BoundaryLocation::YLower};
        case Codim2BoundaryLocation::XUpper_YLower:
            return {BoundaryLocation::XUpper, BoundaryLocation::YLower};
        case Codim2BoundaryLocation::XLower_YUpper:
            return {BoundaryLocation::XLower, BoundaryLocation::YUpper};
        case Codim2BoundaryLocation::XUpper_YUpper:
            return {BoundaryLocation::XUpper, BoundaryLocation::YUpper};

        // X-Z Edges
        case Codim2BoundaryLocation::XLower_ZLower:
            return {BoundaryLocation::XLower, BoundaryLocation::ZLower};
        case Codim2BoundaryLocation::XUpper_ZLower:
            return {BoundaryLocation::XUpper, BoundaryLocation::ZLower};
        case Codim2BoundaryLocation::XLower_ZUpper:
            return {BoundaryLocation::XLower, BoundaryLocation::ZUpper};
        case Codim2BoundaryLocation::XUpper_ZUpper:
            return {BoundaryLocation::XUpper, BoundaryLocation::ZUpper};

        // Y-Z Edges
        case Codim2BoundaryLocation::YLower_ZLower:
            return {BoundaryLocation::YLower, BoundaryLocation::ZLower};
        case Codim2BoundaryLocation::YUpper_ZLower:
            return {BoundaryLocation::YUpper, BoundaryLocation::ZLower};
        case Codim2BoundaryLocation::YLower_ZUpper:
            return {BoundaryLocation::YLower, BoundaryLocation::ZUpper};
        case Codim2BoundaryLocation::YUpper_ZUpper:
            return {BoundaryLocation::YUpper, BoundaryLocation::ZUpper};

        default: throw std::runtime_error("Invalid adjacent boundary location index.");
    }
}

/** @brief Possible locations of a 3-codimensional boundary (a corner in 3D) */
enum class Codim3BoundaryLocation {
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
getAdjacentBoundaryLocations(Codim3BoundaryLocation location)
{
    switch (location)
    {
        // Lower Z Plane
        case Codim3BoundaryLocation::XLower_YLower_ZLower:
            return {BoundaryLocation::XLower, BoundaryLocation::YLower, BoundaryLocation::ZLower};
        case Codim3BoundaryLocation::XUpper_YLower_ZLower:
            return {BoundaryLocation::XUpper, BoundaryLocation::YLower, BoundaryLocation::ZLower};
        case Codim3BoundaryLocation::XLower_YUpper_ZLower:
            return {BoundaryLocation::XLower, BoundaryLocation::YUpper, BoundaryLocation::ZLower};
        case Codim3BoundaryLocation::XUpper_YUpper_ZLower:
            return {BoundaryLocation::XUpper, BoundaryLocation::YUpper, BoundaryLocation::ZLower};

        // Upper Z Plane
        case Codim3BoundaryLocation::XLower_YLower_ZUpper:
            return {BoundaryLocation::XLower, BoundaryLocation::YLower, BoundaryLocation::ZUpper};
        case Codim3BoundaryLocation::XUpper_YLower_ZUpper:
            return {BoundaryLocation::XUpper, BoundaryLocation::YLower, BoundaryLocation::ZUpper};
        case Codim3BoundaryLocation::XLower_YUpper_ZUpper:
            return {BoundaryLocation::XLower, BoundaryLocation::YUpper, BoundaryLocation::ZUpper};
        case Codim3BoundaryLocation::XUpper_YUpper_ZUpper:
            return {BoundaryLocation::XUpper, BoundaryLocation::YUpper, BoundaryLocation::ZUpper};

        default: throw std::runtime_error("Invalid adjacent boundary location index.");
    }
}

/**
 * @brief Get the BoundaryType from input keyword, and throw and error if the keyword does not
 * correspond to any known boundary type.
 */
inline BoundaryType getBoundaryTypeFromString(std::string const& name)
{
    static std::unordered_map<std::string, BoundaryType> const typeMap_ = {
        {"none", BoundaryType::None},
        {"open", BoundaryType::Open},
        {"super-magnetofast-inflow", BoundaryType::SuperMagnetofastInflow},
        {"reflective", BoundaryType::Reflective},
        {"super-magnetofast-outflow", BoundaryType::SuperMagnetofastOutflow},
        {"free-pressure-inflow", BoundaryType::FreePressureInflow},
        {"fixed-pressure-outflow", BoundaryType::FixedPressureOutflow},
    };

    auto it = typeMap_.find(name);
    if (it == typeMap_.end())
        throw std::runtime_error("Wrong boundary type name = " + name);
    return it->second;
}

/**
 * @brief Get the BoundaryType from input keyword, and throw and error if the keyword foes not
 * correspond to any known boundary type.
 */
inline BoundaryLocation getBoundaryLocationFromString(std::string const& name)
{
    static std::unordered_map<std::string, BoundaryLocation> const typeMap_ = {
        {"xlower", BoundaryLocation::XLower}, {"xupper", BoundaryLocation::XUpper},
        {"ylower", BoundaryLocation::YLower}, {"yupper", BoundaryLocation::YUpper},
        {"zlower", BoundaryLocation::ZLower}, {"zupper", BoundaryLocation::ZUpper},
    };

    if (typeMap_.count(name))
        return typeMap_.at(name);
    throw std::runtime_error("Wrong boundary location name = " + name);
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
