#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_DEFS_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_DEFS_HPP

#include <stdexcept>
#include <string>
#include <unordered_map>

namespace PHARE::core
{

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
