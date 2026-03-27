#ifndef PHARE_CORE_INNER_BOUNDARY_SPHERE_INNER_BOUNDARY_HPP
#define PHARE_CORE_INNER_BOUNDARY_SPHERE_INNER_BOUNDARY_HPP

#include <cmath>
#include <stdexcept>
#include <utility>

#include "core/inner_boundary/inner_boundary.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class SphereInnerBoundary : public InnerBoundary<dim>
{
public:
    using Base       = InnerBoundary<dim>;
    using point_type = InnerBoundary<dim>::point_type;

    SphereInnerBoundary(std::string name, point_type center, double radius)
        : Base{name}
        , center_{std::move(center)}
        , radius_{radius}
    {
        if (radius_ <= 0)
            throw std::runtime_error("SphereInnerBoundary radius must be > 0");
    }

    double signedDistance(point_type const& point) const override
    {
        return distance_(point, center_) - radius_;
    }

    point_type normal(point_type const& point) const override
    {
        auto direction = point - center_;
        auto length    = distance_(point, center_);
        if (length == 0.)
            throw std::runtime_error("SphereInnerBoundary normal undefined at center");
        return direction * (1. / length);
    }

private:
    static double distance_(point_type const& a, point_type const& b)
    {
        double squared = 0.;
        for (std::size_t i = 0; i < dim; ++i)
        {
            double const d = a[i] - b[i];
            squared += d * d;
        }
        return std::sqrt(squared);
    }

    point_type center_;
    double radius_;
};

} // namespace PHARE::core

#endif
