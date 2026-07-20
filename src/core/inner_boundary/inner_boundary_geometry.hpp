#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_GEOMETRY_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_GEOMETRY_HPP

#include "core/inner_boundary/inner_boundary_defs.hpp"
#include "core/utilities/point/point.hpp"

#include <limits>

namespace PHARE::core
{

/**
 * @brief Abstract implicit description of an embedded boundary.
 *
 * Concrete implementations provide a signed-distance function and the
 * corresponding outward normal. The default projection and symmetry helpers are
 * derived from that signed-distance representation.
 *
 * @tparam dim Spatial dimension of the embedding space.
 */
template<std::size_t dim>
class InnerBoundaryGeometry
{
public:
    using point_type = Point<double, dim>;

    InnerBoundaryGeometry(std::string name)
        : name_{name}
    {
    }

    virtual ~InnerBoundaryGeometry() = default;

    /**
     * @brief Human-readable boundary name.
     */
    std::string name() const { return name_; }

    /**
     * @brief Boundary shape discriminator.
     */
    InnerBoundaryShape shape() const { return shape_; }

    /**
     * @brief Signed distance to the boundary.
     *
     * Negative values are inside the boundary, positive values are outside, and
     * zero lies on the boundary itself.
     *
     * @param point Physical point where the signed distance is evaluated.
     * @return Signed distance at @p point.
     */
    virtual double signedDistance(point_type const& point) const = 0;

    /**
     * @brief Outward unit normal at the supplied point.
     *
     * @param point Physical point where the normal is evaluated.
     * @return Outward normal direction at @p point.
     */
    virtual point_type normal(point_type const& point) const = 0;

    /**
     * @brief Characteristic length of the boundary, used to decide whether a level
     * resolves it.
     *
     * For a closed convex obstacle this is its smallest radius of curvature (e.g. the
     * sphere radius); for an unbounded boundary (e.g. a plane) it is infinite, since a
     * plane is resolved at any mesh spacing. The inner-boundary classifier compares this
     * to the mesh spacing to switch to a degraded (first-order, ideal) flux near the
     * obstacle on under-resolved levels.
     */
    virtual double characteristicLength() const
    {
        return std::numeric_limits<double>::infinity();
    }

    /**
     * @brief Orthogonal projection of a point onto the boundary.
     *
     * @param point Physical point to project.
     * @return Projection of @p point onto the boundary surface.
     */
    virtual point_type project(point_type const& point) const
    {
        auto n = normal(point);
        return point - n * signedDistance(point);
    }

    /**
     * @brief Mirror image of a point through the boundary surface.
     *
     * @param point Physical point to reflect.
     * @return Symmetric point with respect to the boundary surface.
     */
    virtual point_type symmetric(point_type const& point) const
    {
        auto projected = project(point);
        return projected * 2. - point;
    }

private:
    InnerBoundaryShape shape_;
    std::string name_;
};
} // namespace PHARE::core

#endif
