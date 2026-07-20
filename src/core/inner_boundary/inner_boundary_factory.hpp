#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_FACTORY_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_FACTORY_HPP


#include "core/inner_boundary/inner_boundary_defs.hpp"
#include "core/inner_boundary/inner_boundary_geometry.hpp"
#include "core/inner_boundary/plane_inner_boundary.hpp"
#include "core/inner_boundary/sphere_inner_boundary.hpp"

#include "core/utilities/point/point.hpp"
#include "initializer/data_provider.hpp"

#include <stdexcept>
#include <string>
#include <vector>


namespace PHARE::core
{
template<std::size_t dim>
class InnerBoundaryFactory
{
public:
    using inner_boundary_type = InnerBoundaryGeometry<dim>;
    using point_type          = Point<double, dim>;

    static std::unique_ptr<inner_boundary_type>
    create(initializer::PHAREDict const& simulation_dict)
    {
        if (!simulation_dict.contains("inner_boundary"))
            return nullptr;

        auto const& dict  = simulation_dict["inner_boundary"];
        auto shape_name   = dict["shape"].template to<std::string>();
        auto shape        = getInnerBoundaryShapeFromString(shape_name);
        auto boundaryName = dict["name"].template to<std::string>();

        switch (shape)
        {
            case InnerBoundaryShape::Sphere: {
                auto center = asPoint_(dict["center"].template to<std::vector<double>>());
                auto radius = dict["radius"].template to<double>();
                return std::make_unique<SphereInnerBoundary<dim>>(boundaryName, center, radius);
            }
            case InnerBoundaryShape::Plane: {
                auto point  = asPoint_(dict["point"].template to<std::vector<double>>());
                auto normal = asPoint_(dict["normal"].template to<std::vector<double>>());
                return std::make_unique<PlaneInnerBoundary<dim>>(boundaryName, point, normal);
            }
            default: throw std::runtime_error("Unknow inner boudary shape");
        }
    }

private:
    static point_type asPoint_(std::vector<double> const& coords)
    {
        if (coords.size() != dim)
            throw std::runtime_error("inner_boundary coordinates have invalid dimension");

        point_type point;
        for (std::size_t i = 0; i < dim; ++i)
            point[i] = coords[i];
        return point;
    }
};

} // namespace PHARE::core

#endif
