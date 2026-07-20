#ifndef PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_CONDITION_HPP
#define PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_CONDITION_HPP


#include "core/inner_boundary/inner_boundary_geometry.hpp"

namespace PHARE::core
{

template<size_t dim>
class InnerBoundaryCondition
{
public:
    using inner_boundary_type = InnerBoundaryGeometry<dim>;

    InnerBoundaryCondition() = default;

    virtual ~InnerBoundaryCondition() = default;

private:
};

} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_INNER_BOUNDARY_CONDITION_HPP
