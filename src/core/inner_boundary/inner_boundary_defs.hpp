#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_DEFS_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_DEFS_HPP

#include "core/utilities/point/point.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace PHARE::core
{

/**
 * @brief Status of a mesh element (cell, face, edge, or node) relative to an inner boundary.
 *
 * - **Fluid**    — element lies entirely in the fluid domain.
 * - **Cut**      — element straddles the boundary surface.
 * - **Ghost**    — element lies inside the body but is used to enforce the boundary condition
 *                  by mirror-point interpolation from the fluid side.
 * - **Inactive** — element lies entirely inside the body and plays no role in the solver.
 */
enum class ElemStatus : std::uint8_t { Fluid, Cut, Ghost, Inactive };

/// Convert an ElemStatus value to its double encoding for field storage.
inline constexpr double toDouble(ElemStatus s)
{
    return static_cast<double>(static_cast<std::uint8_t>(s));
}


/**
 * @brief Precomputed per-ghost-element data used by the BC applier every time step.
 *
 * Storing these avoids recomputing expensive boundary queries (normal, symmetric)
 * once per field per ghost element per time step.
 *
 * @note `mirrorIsInterpolable` is set to `false` when the mirror point sits too close
 * to the ghost-box edge for the inner-BC interpolation support (2 consecutive grid
 * values per direction) to fit inside the allocated field extent. When `false`, the
 * BC applier must skip the interpolation.
 *
 * @warning when mirror is not interpolable, nothing is done for the ghost. This might be the cause
 * of issue, TBD. If so, lower order interp could be done.
 * */
template<std::size_t dim>
struct GhostElemData
{
    Point<std::uint32_t, dim> index; ///< Local array index of the ghost element.
    Point<double, dim> mirrorPoint;  ///< Physical coords of the symmetric point in the fluid.
    Point<double, dim> normal;       ///< Unit outward normal at the boundary (ghost → mirror).
    bool mirrorIsInterpolable;       ///< True iff the mirror-point interp supports fits.
};


enum class InnerBoundaryShape { Plane, Sphere };

inline InnerBoundaryShape getInnerBoundaryShapeFromString(std::string const& name)
{
    static std::unordered_map<std::string, InnerBoundaryShape> const shapeMap_{
        {"plane", InnerBoundaryShape::Plane}, {"sphere", InnerBoundaryShape::Sphere}};

    auto it = shapeMap_.find(name);
    if (it == shapeMap_.end())
    {
        throw std::runtime_error("Unknown inner boundary shape " + name);
    }
    return it->second;
}

} // namespace PHARE::core


#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_DEFS_HPP
