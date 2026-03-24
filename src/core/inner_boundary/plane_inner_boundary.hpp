#ifndef PHARE_CORE_INNER_BOUNDARY_PLANE_INNER_BOUNDARY_HPP
#define PHARE_CORE_INNER_BOUNDARY_PLANE_INNER_BOUNDARY_HPP

#include <cmath>
#include <stdexcept>
#include <utility>

#include "core/inner_boundary/inner_boundary.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class PlaneInnerBoundary : public InnerBoundary<dim>
{
public:
    using point_type = typename InnerBoundary<dim>::point_type;

    PlaneInnerBoundary(point_type point, point_type normal)
        : point_{std::move(point)}
        , normal_{normalize_(std::move(normal))}
    {
    }

    double signedDistance(point_type const& point) const override
    {
        return dot_(point - point_, normal_);
    }

    point_type normal(point_type const& /*point*/) const override { return normal_; }

private:
    static double dot_(point_type const& a, point_type const& b)
    {
        double sum = 0.;
        for (std::size_t i = 0; i < dim; ++i)
            sum += a[i] * b[i];
        return sum;
    }

    static point_type normalize_(point_type const& vec)
    {
        double squared = 0.;
        for (std::size_t i = 0; i < dim; ++i)
            squared += vec[i] * vec[i];

        auto const length = std::sqrt(squared);
        if (length == 0.)
            throw std::runtime_error("PlaneInnerBoundary normal cannot be zero");
        return vec * (1. / length);
    }

    point_type point_;
    point_type normal_;
};

} // namespace PHARE::core

#endif
