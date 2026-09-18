#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_DIPOLE_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_DIPOLE_HPP

#include "core/models/external_field_updater_builtin.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>

namespace PHARE::core
{
/**
 * @brief Implements a static dipole external field.
 */
template<typename VecFieldT, typename GridLayoutT>
class ExternalFieldUpdaterDipole
    : public ExternalFieldUpdaterBuiltin<ExternalFieldUpdaterDipole<VecFieldT, GridLayoutT>,
                                         VecFieldT, GridLayoutT>
{
    using Super = ExternalFieldUpdaterBuiltin<ExternalFieldUpdaterDipole<VecFieldT, GridLayoutT>,
                                              VecFieldT, GridLayoutT>;

public:
    using vecfield_type  = Super::vecfield_type;
    using value_type     = Super::value_type;
    using point_type     = Super::point_type;
    using vector_type    = Point<value_type, GridLayoutT::dimension>;
    using component_type = Super::component_type;

    static constexpr std::size_t dimension = Super::dimension;

    /**
     * @param position where the dipole is placed
     * @param moment the moment, one component per dimension
     * @param radius radius of the uniformly magnetized sphere (3D) or cylinder (2D) the dipole
     * is, 0 for a point dipole. Inside it, B0 is constant.
     */
    ExternalFieldUpdaterDipole(point_type position, vector_type moment, double radius)
        : position_{position}
        , moment_{moment}
        , radiusSquared_{radius * radius}
    {
        if (!std::isfinite(radius) or radius < 0. or !std::isfinite(radiusSquared_))
            throw std::invalid_argument("dipole radius must be finite and positive or zero, got "
                                        + std::to_string(radius));
    }

    template<component_type i>
    double potential(point_type const& coords, double /*time*/) const
    {
        if constexpr (dimension == 1)
        {
            throw std::runtime_error("a dipole in 1D makes no sense");
        }
        else if constexpr (dimension == 2)
        {
            if constexpr (i == component_type::Z)
            {
                constexpr double factor = 1. / (2. * std::numbers::pi);
                point_type const r      = coords - position_;
                double const rSquared   = std::max(
                    std::inner_product(r.begin(), r.end(), r.begin(), 0.0), radiusSquared_);
                if (rSquared == 0.)
                    return 0.0;
                // z component of the cross product `moment_` times `r`
                return factor * (moment_[0] * r[1] - moment_[1] * r[0]) / rSquared;
            }
            else
                return 0.0;
        }
        else // 3D case
        {
            point_type const r = coords - position_;
            double const rSquared
                = std::max(std::inner_product(r.begin(), r.end(), r.begin(), 0.0), radiusSquared_);
            if (rSquared == 0.)
                return 0.0;
            double constexpr factor = 1. / (4. * std::numbers::pi);
            // elegant trick to express component i of cross product `moment_` times `r`
            constexpr auto j = (static_cast<std::size_t>(i) + 1) % 3;
            constexpr auto k = (static_cast<std::size_t>(i) + 2) % 3;
            return factor * (moment_[j] * r[k] - moment_[k] * r[j])
                   / (rSquared * std::sqrt(rSquared));
        }
    }

private:
    point_type position_;  //!< where is placed the dipole in space
    vector_type moment_;   //!< the moment vector, one component per dimension
    double radiusSquared_; //!< below this squared distance to the dipole, B0 is uniform
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_DIPOLE_HPP
